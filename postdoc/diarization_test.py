import os
from pydub import AudioSegment
from pyAudioAnalysis import audioSegmentation as aS
from tqdm import tqdm
import time
import subprocess

def extract_audio(video_path, audio_path, n_seconds):
    command = f"ffmpeg -i {video_path} -vn -acodec pcm_s16le -ar 44100 -ac 2 -t {n_seconds} {audio_path}"
    subprocess.call(command, shell=True)

def diarize_audio(file_path):
    [flags, classes, acc, segs] = aS.speaker_diarization(file_path, 4, plot_res=False)
    return segs, classes

def extract_segments(audio_path, segments, classes):
    audio_full = AudioSegment.from_file(audio_path)
    extracted_segments = []
    base_dir = os.path.dirname(audio_path)
    segments_dir = os.path.join(base_dir, "segments")
    os.makedirs(segments_dir, exist_ok=True)  # Create directory for segments
    
    manifest_path = os.path.join(base_dir, "segments_manifest.txt")
    with open(manifest_path, 'w') as manifest_file:
        for i, (start, end) in enumerate(segments):
            segment = audio_full[start*1000:end*1000]  # pydub uses milliseconds
            segment_file_name = f"segment_{i}_class_{classes[i]}.wav"
            segment_file_path = os.path.join(segments_dir, segment_file_name)
            segment.export(segment_file_path, format="wav")
            extracted_segments.append(segment_file_path)

            # Write segment info to manifest
            manifest_file.write(f"{segment_file_name}\t{start}\t{end}\tClass: {classes[i]}\n")

    return extracted_segments

def main():

    # general purpose variables
    video_path = "~/test.mp4"
    video_path = os.path.expanduser(video_path)  # Expands the '~' to the full home directory path
    audio_path = "~/test_audio.wav"
    audio_path = os.path.expanduser(audio_path)  # Expands the '~' to the full home directory path

    #get audio from video
    print("Starting audio extraction")
    extract_audio(video_path, audio_path, 150)

    #perform diarization
    print("Starting diarization...")
    start_time = time.time()
    segments, classes = diarize_audio(audio_path)
    end_time = time.time()
    print(f"Diarization completed in {end_time - start_time:.2f} seconds.")

    #extract and save segments
    print("Extracting and saving segments...")
    for segment in tqdm(extract_segments(audio_path, segments, classes), desc="Processing segments"):
        pass  # The actual work is done inside the extract_segments function

    print("All segments processed and saved.")

if __name__ == "__main__":
    main()
