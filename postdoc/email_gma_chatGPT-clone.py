import smtplib
import time
import imaplib
import email
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart  
from email.mime.image import MIMEImage
from email import policy  # Import policy here
from openai import OpenAI
import json
import re
import requests
import os


# Initialize OpenAI client globally
client = OpenAI()

# Email credentials and settings
imap_server = "imap.gmail.com"
smtp_server = "smtp.gmail.com"
email_address = os.environ.get('EMAIL_APP')
password = os.environ.get('EMAIL_APP_PASSWORD')

def save_conversations():
    with open('conversation_histories.json', 'w') as file:
        json.dump(conversation_histories, file)

def load_conversations():
    try:
        with open('conversation_histories.json', 'r') as file:
            return json.load(file)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}  # Return an empty dictionary if the file does not exist or is empty

conversation_histories = load_conversations()

import os

def check_inbox():
    mail = imaplib.IMAP4_SSL(imap_server)
    mail.login(email_address, password)
    mail.select('inbox')

    status, data = mail.search(None, 'UNSEEN')
    mail_ids = []

    if status == 'OK':
        mail_ids = data[0].split()
    
    if mail_ids:
        latest_email_id = mail_ids[-1]
        status, data = mail.fetch(latest_email_id, '(RFC822)')
        raw_email = data[0][1]
        email_message = email.message_from_bytes(raw_email, policy=policy.default)
        sender_email = email.utils.parseaddr(email_message['From'])[1]

        message_body = ""
        image_data = None
        if email_message.is_multipart():
            for part in email_message.walk():
                content_type = part.get_content_type()
                content_disposition = str(part.get("Content-Disposition"))
                
                if "attachment" in content_disposition:
                    filename = part.get_filename()
                    if filename and filename.endswith(('.png', '.jpg', '.jpeg', '.gif')):
                        # Save image data; you can also save the file if needed
                        image_data = part.get_payload(decode=True)
                        print(f"Image found in attachment: {filename}")
                        # Optionally, save the image file
                        with open(os.path.join('images', filename), 'wb') as f:
                            f.write(image_data)

                    elif filename == "text_0.txt":
                        message_body = part.get_payload(decode=True).decode('utf-8')
                        print("Message found in attachment.")
                        break

                elif content_type == "text/plain" and "attachment" not in content_disposition:
                    message_body = part.get_payload(decode=True).decode('utf-8')
                    print("Message found in email body.")
                    break

        if message_body:
            return message_body, sender_email, image_data
        else:
            print("No message body found.")
            return None, None, None
    else:
        return None, None, None


def process_message(message_text, sender_email, image_data=None, image_filename=None):

    # Send the response text back to the sender
    grandma_phone_number = "6028199984"
    immediate_response_english = "Got your message! Let me think for a sec."
    immediate_response_russian = "Получил ваше сообщение! Дайте мне секундочку подумать."

    # Check if the message is from grandma's phone number
    if sender_email.endswith(grandma_phone_number + "@vzwpix.com"):
        # Send an immediate response in Russian
        send_response(immediate_response_russian, sender_email)
    else:
        # Send an immediate response in English
        send_response(immediate_response_english, sender_email)
        
    # Custom instructions for AI to suggest images and keep things brief
    custom_instructions = (
    "If an image would help explain or enhance the response, include the tag '[[Generate Image]]' in double square brackets followed immediately by a concise, yet detailed image description contained in curly braces, like so: '{...}'. Commas may separate different elements of the image description. Never include more than a single [[Generate Image]] tag with paired image description in {...} per response. As an AI assistant, you may not be able to generate images directly, but you can include these tags for an independent program to generate a requested image. Please put the image information at the very end of your response. "
    "Use relevant examples for clarification. "
    "Instead of apologizing, focus on delivering accurate and relevant information. "
    "If events or information are beyond your scope or knowledge cutoff date in September 2021, provide a response stating 'I don't know'. "
    "Refrain from disclaimers about you not being a professional or expert. "
    "Keep responses unique and free of repetition. "
    "Break down complex problems or tasks into smaller, manageable steps and think through each one, sharing your reasoning as you go. "
    "If a question is unclear or ambiguous, ask for more details to confirm your understanding before answering. "
    "Aim to be concise and prioritize making every text response count, omitting needless words. "
    )

    # Custom instructions for grandma
    grandma_profile = "Вы общаетесь с пожилым человеком, родившимся в декабре 1937 года, который говорит только по-русски и живет в Седоне, штат Аризона. Она была неврологом и может задавать вопросы, связанные с технической поддержкой или другими темами. Важно четко и просто отвечать на ее вопросы на русском языке. Если для ответа не хватает информации, пожалуйста, задавайте уточняющие вопросы, чтобы помочь ей предоставить необходимые детали."


    # Check if the message is from grandma's phone number
    if sender_email.endswith(grandma_phone_number + "@vzwpix.com"):
        profile_context = "Profile Context: " + grandma_profile
    else:
        profile_context = ""

    # Retrieve the specific conversation history for this sender
    conversation = conversation_histories.get(sender_email, [])

    # Add the received message to the conversation history
    conversation.append({"role": "user", "content": message_text})

    # Add an image attachment line if image data is present
    # TODO implement rehosting of image on aws or wherever to pass to gpt-4-vision-preview
    # if image_data:
    interpret_image = False
    if interpret_image:
    
        # Add a reference to the image in the conversation
        conversation.append({"role": "user", "content": f"attachment: {image_filename}"})
        
        # Use GPT-4 Vision model for image processing
        image_model = "gpt-4-vision-preview"
        image_messages = [
            {
                "role": "user",
                "content": [
                    {"type": "text", "text": message_text},
                    {"type": "image_url", "image_url": {"url": image_data}}
                ],
            }
        ]

        # Call to OpenAI API for image processing
        image_completion = client.chat.completions.create(
            model=image_model,
            messages=image_messages
        )

        # Extract the description from the image_completion response
        image_description = image_completion.choices[0].message.content if hasattr(image_completion.choices[0].message, 'content') else str(image_completion.choices[0].message)

        # Append the image description to the conversation history
        conversation.append({"role": "user", "content": f"User image description: {image_description}"})
        
    # Specify model eg GPT-4 or GPT-3.5-Turbo for text processing
    # model = "gpt-3.5-turbo"
    model = "gpt-4"
        
    # Construct the message for text processing
    truncated_conversation = truncate_to_recent_tokens(conversation)
    messages = [{"role": "system", "content": custom_instructions + "\n\n" + profile_context}] + truncated_conversation

    # Initialize OpenAI client and send the conversation history
    completion = client.chat.completions.create(
        model=model,
        messages=messages
    )

    # Get the response from OpenAI
    response_object = completion.choices[0].message
    response_text = response_object.content if hasattr(response_object, 'content') else str(response_object)

    # Check for image generation marker
    image_url = None
    if "[[Generate Image]]" in response_text:
        image_description = extract_image_description(response_text)
        if image_description:
            image_url = generate_image_with_dall_e(image_description)

    # Add the AI's response to the conversation history
    conversation.append({"role": "assistant", "content": response_text})
    conversation_histories[sender_email] = conversation

    # Save the updated conversations
    save_conversations()

    print("Response message: ", response_text)

    # Remove image caption
    response_text = remove_image_marker(response_text)

    # Send the response text back to the sender
    send_response(response_text, sender_email)
    
    if image_url:
        send_image_response(image_url, sender_email)

def truncate_to_recent_tokens(conversation, max_tokens=2048):
    total_tokens = 0
    truncated_conversation = []

    # Iterate over the conversation in reverse (starting from the most recent)
    for message in reversed(conversation):
        message_content = message["content"]
        message_tokens = len(message_content.split())  # Simple tokenization by whitespace

        if total_tokens + message_tokens > max_tokens:
            break  # Stop if adding this message would exceed the token limit

        truncated_conversation.append(message)
        total_tokens += message_tokens

    return list(reversed(truncated_conversation))  # Reverse again to maintain chronological order


# Function to generate an image using DALL-E
def generate_image_with_dall_e(image_description):
    response = client.images.generate(
        model="dall-e-3",
        prompt=image_description,
        size="1024x1024",
        quality="standard",
        n=1,
    )
    if response.data:
        return response.data[0].url
    else:
        return None


def send_image_response(image_url, recipient_email):
    sender_email = "robot.kolka.vetr@gmail.com"
    
    # Create MIME multipart message
    msg = MIMEMultipart()
    msg['From'] = sender_email
    msg['To'] = recipient_email
    msg['Subject'] = ""

    # Download the image from the URL
    response = requests.get(image_url)
    if response.status_code == 200:
        image_data = response.content
        image = MIMEImage(image_data)
        msg.attach(image)
    else:
        print("Failed to download image.")
        return

    # Connect to Gmail's SMTP server
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(sender_email, password)

    # Send the email
    server.send_message(msg)
    server.quit()

    print("Image response sent to", recipient_email)


# Function to extract image description and remove the marker
def extract_image_description(response_text):
    # Pattern to match the image description structure {~{...}~}
    pattern = r"\{(.+?)\}"
    match = re.search(pattern, response_text)
    
    if match:
        # Extract the image description
        return match.group(1)
    else:
        # Handle the case where the image description is not found
        return None

def remove_image_marker(response_text):
    # Pattern to match and remove the entire structure including the marker
    pattern = r"\[\[Generate Image]\].*?\{.+?\}"
    modified_response_text = re.sub(pattern, '', response_text).strip()

    return modified_response_text

def send_response(response, recipient_email):

    # Sender's email and password (use app password if 2FA is enabled)
    sender_email = "robot.kolka.vetr@gmail.com"

    # Create MIME multipart message
    msg = MIMEMultipart()
    msg['From'] = sender_email  # The bot's email address
    msg['To'] = recipient_email  # The recipient's email address
    msg['Subject'] = "Робо-Колька тебе отвечает:"

    # Attach the response text
    msg.attach(MIMEText(response, 'plain'))

    # Connect to Gmail's SMTP server
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()  # Secure the connection
    server.login(sender_email, password)

    # Send the email
    server.send_message(msg)
    server.quit()

    print("Response sent to", recipient_email)



def main():
    print("Starting the script...")
    while True:
        message, sender_email, image = check_inbox()
        if message:
            print(f"Processing message: {message}")
            process_message(message, sender_email)
        time.sleep(10)
            
if __name__ == "__main__":
    main()
