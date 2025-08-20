import os
import re
from glob import glob
from datetime import datetime, timedelta
import numpy as np
from openai import OpenAI

def read_files(directory_path):
    files = sorted(glob(os.path.join(directory_path, "*_*.txt")), key=os.path.getmtime)
    transcripts = []

    for file in files:
        with open(file, 'r') as f:
            transcripts.append(f.read())

    return transcripts
def parse_transcripts(transcripts):
    speech_pattern = re.compile(r"(\w+)\s+(\d+:\d+(:\d+)?)(.*)")
    parsed_data = []

    current_speaker = None
    current_timestamp = None
    current_text = []

    for transcript in transcripts:
        for line in transcript.split('\n'):
            match = speech_pattern.match(line)
            if match:
                # If there is a current speaker, save their accumulated text
                if current_speaker and current_timestamp:
                    parsed_data.append((current_speaker, current_timestamp, ' '.join(current_text).strip()))

                # Start a new entry
                current_speaker = match.group(1)
                current_timestamp = match.group(2)
                current_text = [match.group(4).strip()] if match.group(4) else []

            else:
                if current_text is not None:
                    current_text.append(line.strip())

    # Don't forget to save the last speaker's accumulated text
    if current_speaker and current_timestamp:
        parsed_data.append((current_speaker, current_timestamp, ' '.join(current_text).strip()))

    return parsed_data

def sort_by_timestamp(parsed_data):
    def convert_timestamp(ts):
        parts = ts.split(':')
        parts = list(map(int, parts))
        return datetime.strptime(':'.join(map(str, parts)), '%H:%M:%S' if len(parts) == 3 else '%M:%S')

    sorted_data = sorted(parsed_data, key=lambda x: convert_timestamp(x[1]))
    return sorted_data

def adjust_timestamps(parsed_data):
    def convert_timestamp(ts):
        parts = ts.split(':')
        parts = list(map(int, parts))
        return timedelta(hours=parts[0], minutes=parts[1], seconds=parts[2] if len(parts) == 3 else 0)

    def timestamp_to_str(td):
        total_seconds = int(td.total_seconds())
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if hours > 0:
            return f"{hours}:{minutes:02}:{seconds:02}"
        else:
            return f"{minutes}:{seconds:02}"

    adjusted_data = []
    cumulative_offset = timedelta()
    previous_time = None

    for entry in parsed_data:
        speaker, timestamp, text = entry
        current_time = convert_timestamp(timestamp)

        if previous_time is not None and current_time < previous_time:
            cumulative_offset += timedelta(minutes=1)

        adjusted_time = current_time + cumulative_offset
        adjusted_timestamp = timestamp_to_str(adjusted_time)
        adjusted_data.append((speaker, adjusted_timestamp, text))
        previous_time = adjusted_time

    return adjusted_data

def format_output(sorted_data):
    formatted_output = []
    current_speaker = None

    for entry in sorted_data:
        speaker, timestamp, text = entry
        if speaker != current_speaker:
            formatted_output.append(f"\n{speaker} {timestamp}\n{text}")
            current_speaker = speaker
        else:
            formatted_output[-1] += f" {text}"

    return '\n'.join(formatted_output)

def adjust_timestamps(parsed_data):
    # Step 1: Extract the timestamps
    timestamps = [entry[1] for entry in parsed_data]

    # Step 2: Convert timestamps to timedelta objects
    def convert_timestamp(ts):
        parts = list(map(int, ts.split(':')))
        if len(parts) == 2:
            return timedelta(minutes=parts[0], seconds=parts[1])
        elif len(parts) == 3:
            return timedelta(hours=parts[0], minutes=parts[1], seconds=parts[2])
        else:
            raise ValueError("Invalid timestamp format")

    timedeltas = [convert_timestamp(ts) for ts in timestamps]

    # Step 3: Find the differences between adjacent timestamps
    differences = [(timedeltas[i+1] - timedeltas[i]).total_seconds() for i in range(len(timedeltas)-1)]

    # Step 4: Set all negative differences to 60 seconds (1 minute)
    differences = [diff if diff >= 0 else 60 for diff in differences]

    # Step 5: Compute the cumulative sum of these differences, appending 0 to the start
    cumulative_sum = np.cumsum([0] + differences)

    # Step 6: Convert the cumulative sums back to timestamps
    def seconds_to_timestamp(seconds):
        total_seconds = int(seconds)
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if hours > 0:
            return f"{hours}:{minutes:02}:{seconds:02}"
        else:
            return f"{minutes}:{seconds:02}"

    new_timestamps = [seconds_to_timestamp(seconds) for seconds in cumulative_sum]

    # Step 7: Reassign this new vector to the parsed_data object in the timestamp field
    adjusted_data = [(entry[0], new_timestamps[i], entry[2]) for i, entry in enumerate(parsed_data)]

    return adjusted_data


def save_chunks(adjusted_data, output_dir, chunk_length=600):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    current_chunk = []
    current_duration = timedelta()
    chunk_index = 1
    previous_entry_time = timedelta()

    # Function to convert timestamp to timedelta
    def convert_timestamp(ts):
        parts = list(map(int, ts.split(':')))
        if len(parts) == 2:
            return timedelta(minutes=parts[0], seconds=parts[1])
        elif len(parts) == 3:
            return timedelta(hours=parts[0], minutes=parts[1], seconds=parts[2])
        else:
            raise ValueError("Invalid timestamp format")

    for entry in adjusted_data:
        speaker, timestamp, text = entry
        entry_time = convert_timestamp(timestamp)

        # Accumulate the current duration with the difference from the previous entry
        if current_duration + (entry_time - previous_entry_time) > timedelta(seconds=chunk_length):
            # Save the current chunk to a file
            chunk_filename = os.path.join(output_dir, f"chunk_{chunk_index}.txt")
            with open(chunk_filename, 'w') as chunk_file:
                for chunk_entry in current_chunk:
                    chunk_file.write(f"{chunk_entry[0]} {chunk_entry[1]}\n{chunk_entry[2]}\n\n")
            chunk_index += 1

            # Start a new chunk with the current entry
            current_chunk = [entry]
            current_duration = entry_time - previous_entry_time
        else:
            current_chunk.append(entry)
            current_duration += entry_time - previous_entry_time

        previous_entry_time = entry_time

    # Save the last chunk if it has any entries
    if current_chunk:
        chunk_filename = os.path.join(output_dir, f"chunk_{chunk_index}.txt")
        with open(chunk_filename, 'w') as chunk_file:
            for chunk_entry in current_chunk:
                chunk_file.write(f"{chunk_entry[0]} {chunk_entry[1]}\n{chunk_entry[2]}\n\n")

def initialize_openai_api():
    api_key = os.getenv("OPENAI_API_KEY")
    client = OpenAI(api_key=api_key)
    return client

def generate_meeting_minutes(client, transcript_chunk, previous_minutes, instructions, model="gpt-4o-mini"):
    messages = [
        {"role": "system", "content": instructions},
        {"role": "user", "content": f"Previous minutes:\n{previous_minutes}\n\nCurrent transcript chunk:\n{transcript_chunk}"}
    ]

    response = client.chat.completions.create(
        model=model,
        messages=messages
    )

    return response.choices[0].message.content

def process_chunks_to_minutes(chunks_directory, output_directory, chunk_length=600):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Initialize OpenAI API client
    client = initialize_openai_api()

    # Define the instructions
    general_instructions = (
        "You are an AI assistant whose job is to convert plaintext transcript chunks into formal meeting minutes that are nicely "
        "formatted in markdown. Ensure that each section is clear and concise, with proper headings and bullet points in simple, straightforward language. "
        "If multiple topics are discussed under a similar theme, "
        "group them under subheadings and sub-subheadings. Follow the minutes format strictly and ensure continuity across chunks. To help you, we will provide the "
        "previous set of meeting minutes, as well as the current transcript chunk. You should be detailed, where appropriate, but  "
        "strive to broadly summarize the contents of each transcript chunk in general terms, instead of describing "
        "everything that was said, and sparingly include details, only when they are very important. Adhere to best practice for writing minutes."
        "Unless otherwise noted, when generating minutes from each chunk, please continue on seamlessly from the previous set of minutes, and end each set of minutes in anticipation of seamless continuation from the previous section."
        "Also, avoid sharing any sensitive information from the minutes. If specific employees' difficulties or affairs are mentioned, omit those details, or summarize them only at a very high level. "
        "In your response, please provide only meeting minutes for the current transcript chunk as a direct continuation of the previous chunk, and also anticipating direct continuation by the following chunk, as they will all be concatenated into a single document after completion."
    )

    first_chunk_instructions = (
        "This is the first transcript chunk, so you may start a new minutes documentfor the quarterly meeting of the Board of Directors as the Wild Animal Initiative."
    )

    # Get the list of chunk files
    chunk_files = sorted(glob(os.path.join(chunks_directory, "chunk_*.txt")), key=os.path.getmtime)
    previous_minutes = ""

    for i, chunk_file in enumerate(chunk_files):
        with open(chunk_file, 'r') as file:
            transcript_chunk = file.read()

        if i == 0:
            instructions = first_chunk_instructions + "\n\n" + general_instructions
        else:
            instructions = general_instructions

        minutes = generate_meeting_minutes(client, transcript_chunk, previous_minutes, instructions)

        output_file = os.path.join(output_directory, f"minutes_chunk_{i+1}.md")
        with open(output_file, 'w') as out_file:
            out_file.write(minutes)

        previous_minutes += minutes + "\n\n"

def concatenate_minutes_files(minutes_directory, output_file):
    # Get the list of minutes files
    minutes_files = sorted(glob(os.path.join(minutes_directory, "minutes_chunk_*.md")), key=os.path.getmtime)
    
    with open(output_file, 'w') as outfile:
        for minutes_file in minutes_files:
            with open(minutes_file, 'r') as infile:
                outfile.write(infile.read())
                outfile.write("\n\n")  # Add spacing between chunks

def convert_markdown_to_html(markdown_file, output_file):
    subprocess.run(['pandoc', markdown_file, '-o', output_file], check=True)
    print(f"Converted {markdown_file} to {output_file}")


def improve_minutes(concatenated_minutes_file, improved_minutes_file, instructions, model="gpt-4o-mini"):
    # Initialize OpenAI API client
    client = initialize_openai_api()

    # Read the content of the concatenated minutes file
    with open(concatenated_minutes_file, 'r') as file:
        concatenated_minutes = file.read()

    # Define the improvement instructions
    improvement_instructions = (
        instructions + "\n\n"
        "A previous GPT has generated a set of meeting minutes, but it was very wordy in its output, and could not refrain from including things like (continuation) in unnecessary places. Please improve the minutes it provided, rewriting them to be more concise, removing any unnecessary verbosity, and improving their overall flow. Make sure to remove any redundancy between sections, and format your output as markdown with sections and subsections. Please ensure the resulting output is between 3-5 pages long -- so not too brief, but not too wordy, in easy-to-read, understandable, natural language, and organized with bullet points."
    )

    messages = [
        {"role": "system", "content": improvement_instructions},
        {"role": "user", "content": f"Here are the meeting minutes that need improvement:\n\n{concatenated_minutes}"}
    ]

    response = client.chat.completions.create(
        model=model,
        messages=messages
    )

    improved_minutes = response.choices[0].message.content

    # Save the improved minutes to a new file
    with open(improved_minutes_file, 'w') as file:
        file.write(improved_minutes)

    print(f"Improved minutes saved to {improved_minutes_file}")


def main():
    directory_path = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04")
    transcripts = read_files(directory_path)
    parsed_data = parse_transcripts(transcripts)
    adjusted_data = adjust_timestamps(parsed_data)
    chunks_directory = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04/chunks")
    save_chunks(adjusted_data, chunks_directory, 1200)
    minutes_directory = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04/minutes")
    process_chunks_to_minutes(chunks_directory, minutes_directory)
    
    concatenated_minutes_file = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04/complete_minutes.md")
    improved_minutes_file = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04/improved_minutes.md")
    minutes_file_html = os.path.expanduser("~/Documents/Documents - nikolai/WAI_minutes/2024-04/complete_minutes.html")
    concatenate_minutes_files(minutes_directory, concatenated_minutes_file)
    improve_minutes(concatenated_minutes_file, improved_minutes_file, "")
    convert_markdown_to_html(improved_minutes_file, minutes_file_html)

if __name__ == "__main__":
    main()


