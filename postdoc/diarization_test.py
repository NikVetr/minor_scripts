import os
from pydub import AudioSegment
from pyAudioAnalysis import audioSegmentation as aS
from tqdm import tqdm
import time
import subprocess
import shutil
import numpy as np
from scipy.stats import norm
from collections import defaultdict

def preprocess_audio(input_path, output_path):
    """
    Preprocesses an audio file by applying volume normalization and converting it to a suitable format.
    Adjusts the silence removal parameters and sets ffmpeg loglevel to 'error'.
    Automatically overwrites existing files without prompt.
    """
    # Step 1: Volume normalization with automatic overwrite
    normalized_audio = output_path.replace('.wav', '_normalized.wav')
    subprocess.run([
        'ffmpeg', '-y', '-loglevel', 'error', '-i', input_path, '-filter:a', 'loudnorm', normalized_audio
    ], check=True)

    # Step 3: Convert to WAV format, 16000 Hz, mono, with automatic overwrite
    subprocess.run([
        'ffmpeg', '-y', '-loglevel', 'error', '-i', normalized_audio, '-ar', '16000', '-ac', '1', output_path
    ], check=True)

    # Clean up intermediate files
    subprocess.run(['rm', normalized_audio], check=True)
    print(f"Preprocessed audio saved to {output_path}")

def extract_audio(video_path, audio_path, n_seconds):
    command = f"ffmpeg -y -loglevel error -i {video_path} -vn -acodec pcm_s16le -ar 44100 -ac 2 -t {n_seconds} {audio_path}"
    subprocess.call(command, shell=True)

def diarize_audio(file_path, num_speakers=0):
    labels, purity_cluster_m, purity_speaker_m = aS.speaker_diarization(file_path, num_speakers, lda_dim=0, plot_res=False)
    if purity_cluster_m < 0.9:
        print(f"Diarization cluster purity is low, purity_cluster_m = {purity_cluster_m}")
    if purity_speaker_m < 0.9:
        print(f"Diarization speaker purity is low, purity_speaker_m = {purity_speaker_m}")
        
    return labels, purity_cluster_m, purity_speaker_m

def merge_short_segments(labels, window_width=600, weights_type="gaussian"):
    if weights_type == "triangular":
        half = window_width // 2
        weights = np.linspace(1, half, half) / half
        weights = np.concatenate((weights[::-1], [1], weights))  # Symmetric weights with peak at center
    elif weights_type == "gaussian":
        x = np.linspace(-norm.ppf(0.999), norm.ppf(0.999), window_width // 2)
        weights = norm.pdf(x)
        weights = np.concatenate((weights[::-1], [1], weights))  # Ensure symmetry
        weights /= weights.max()  # Normalize so center weight is 1

    labels_np = np.array(labels, dtype=int)
    merged_labels = np.zeros_like(labels_np)

    for i in range(len(labels_np)):
        start_idx = max(i - window_width // 2, 0)
        end_idx = min(i + window_width // 2 + 1, len(labels_np))

        # Adjust weights for the current window
        current_weights = weights[window_width // 2 - i + start_idx : window_width // 2 + end_idx - i]
        window_labels = labels_np[start_idx:end_idx]

        # Compute weighted counts for each label in the window
        unique_labels, counts = np.unique(window_labels, return_counts=True)
        label_weights = np.array([np.sum(current_weights[window_labels == label]) for label in unique_labels])

        # Choose the label with the maximum weighted sum
        max_weight_label = unique_labels[np.argmax(label_weights)]
        merged_labels[i] = max_weight_label

    return merged_labels.tolist()

def extract_segments(audio_path, labels, purity_cluster_m, min_length = 20):
    frame_duration = 0.1  # Duration of each frame in seconds, where each label represents 0.1s
    audio_full = AudioSegment.from_file(audio_path)
    extracted_segments = []
    base_dir = os.path.dirname(audio_path)
    segments_dir = os.path.join(base_dir, "segments")

    # Check if the segments directory exists and clear it
    if os.path.exists(segments_dir):
        shutil.rmtree(segments_dir)
    os.makedirs(segments_dir, exist_ok=True)

    manifest_path = os.path.join(base_dir, "segments_manifest.txt")
    with open(manifest_path, 'w') as manifest_file:
        start_frame = 0
        current_speaker = labels[0]

        for i, speaker in enumerate(labels):
            if speaker != current_speaker or i == len(labels) - 1:
                end_frame = i
                segment_duration = (end_frame - start_frame)
                
                # Check if the segment meets the minimum length requirement
                if segment_duration >= min_length:
                    start_time = start_frame * frame_duration * 1000  # Convert to milliseconds
                    end_time = end_frame * frame_duration * 1000
                    segment_file_name = f"segment_{start_frame}_{end_frame}_class_{current_speaker}.wav"
                    segment_file_path = os.path.join(segments_dir, segment_file_name)
                    segment = audio_full[start_time:end_time]
                    segment.export(segment_file_path, format="wav")
                    extracted_segments.append(segment_file_path)
                    manifest_file.write(f"{segment_file_name}\t{start_time/1000}\t{end_time/1000}\tClass: {current_speaker}\n")
                
                start_frame = i
                current_speaker = speaker

    return extracted_segments

def main():

    # general purpose variables
    video_path = "~/test.mp4"
    video_path = os.path.expanduser(video_path)  # Expands the '~' to the full home directory path
    raw_audio_path = "~/raw_test_audio.wav"
    raw_audio_path = os.path.expanduser(raw_audio_path)  # Expands the '~' to the full home directory path
    audio_path = "~/test_audio.wav"
    audio_path = os.path.expanduser(audio_path)  # Expands the '~' to the full home directory path

    #get audio from video
    print("Starting audio extraction.")
    extract_audio(video_path, raw_audio_path, 25000)
    
    #clean audio for diarization
    print("Cleaning audio file.")
    preprocess_audio(raw_audio_path, audio_path)

    #perform diarization
    print("Starting diarization...")
    start_time = time.time()
    raw_labels, purity_cluster_m, purity_speaker_m = diarize_audio(audio_path, 3)
    end_time = time.time()
    print(f"Diarization completed in {end_time - start_time:.2f} seconds.")

    # Merge short segments
    labels = merge_short_segments(raw_labels, 500)
    
    #extract and save segments
    print("Extracting and saving segments...")
    for segment in tqdm(extract_segments(audio_path, labels, purity_cluster_m, min_length = 20), desc="Processing segments"):
        pass  # The actual work is done inside the extract_segments function

    print("All segments processed and saved.")

if __name__ == "__main__":
    main()




