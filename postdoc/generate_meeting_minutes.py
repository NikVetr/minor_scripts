import subprocess
import os
from whisper.transcribe import transcribe
import openai
from whisper import load_model

def transcribe_audio(audio_path, model_size="medium"):
    # Load the Whisper model
    model = load_model(model_size)
    # Transcribe the audio
    result = model.transcribe(audio_path)
    return result["text"]

def extract_audio(video_path, audio_path):
    command = f"ffmpeg -i {video_path} -vn -acodec libmp3lame -q:a 2 {audio_path}"

    subprocess.call(command, shell=True)

def summarize_transcript(transcript):
    openai.api_key = os.getenv("OPENAI_API_KEY")
    response = openai.Completion.create(
        model="text-davinci-003",
        prompt=f"Please produce appropriately formatted meeting minutes for the following meeting transcript:\n\n{transcript}",
        max_tokens=1024
    )
    return response.choices[0].text.strip()

def main():
    video_path = "~/test.mp4"
    video_path = os.path.expanduser(video_path)  # Expands the '~' to the full home directory path

    audio_path = "~/test_audio.aac"
    audio_path = os.path.expanduser(audio_path)  # Expands the '~' to the full home directory path

    extract_audio(video_path, audio_path)
    transcript = transcribe_audio(audio_path)
    print("Meeting Summary:")
    print(transcript)

    #summary = summarize_transcript(transcript)
    #print("Meeting Summary:")
    #print(summary)

if __name__ == "__main__":
    main()

