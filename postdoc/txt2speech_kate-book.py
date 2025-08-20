#!/usr/bin/env python3
# tts_batch_parallel.py  — run inside Positron

from pathlib import Path
import io, time, shutil, textwrap, concurrent.futures as cf
from typing import Iterable, Sequence

import tiktoken
import openai
from openai import OpenAI, APIConnectionError
from httpx   import RemoteProtocolError, ReadError
from pydub   import AudioSegment

# ─────────────────────────────  style  ──────────────────────────────
MODEL  = "gpt-4o-mini-tts"
VOICE  = "shimmer"
INSTRUCTIONS = textwrap.dedent("""\
    Voice affect: soft, enveloping, unintrusive—like a gentle ambient presence, softening all t sounds.
    Tone: consistently calm with a narrow dynamic range; avoid sudden inflection.
    Pacing: slow and steady, each word drifting into the next for a flowing cadence.
    Emotion: warm, deeply soothing, subtly nurturing—evoke quiet reassurance.
    Pronunciation: smooth and rounded; soften consonants, slightly elongate vowels.
    Pauses: unhurried and thoughtfully placed, never breaking the serene flow.
    Volume: Consistently low, gently modulated to avoid any startling shifts, maintaining a comforting intimacy.
    Energy: Soft and mellow, with a calm, subdued energy, no engagement or excitement.
    Cadence: Slow and fluid, with a naturally rolling rhythm, one word flows into the next""")
# ────────────────────────────────────────────────────────────────────

# ───────────────────────────  folders  ──────────────────────────────
#SRC_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/the_keys/sections")
#SRC_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_creatures_great_and_small/sections/")
#SRC_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/nothing_much_happens/sections/")
SRC_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_things_wise_and_wonderful/sections/")
#SRC_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_things_bright_and_beautiful/sections/")
OUT_DIR = SRC_DIR / VOICE / MODEL / "mp3"
# ────────────────────────────────────────────────────────────────────

# ─────────────  token-aware text chunker  ─────────────
ENC = tiktoken.encoding_for_model(MODEL)
TOK_LIM = 1_800                       # leave margin under 2 000

def tokens(txt: str) -> int:
    return len(ENC.encode(txt))

def chunk_text(txt: str) -> list[str]:
    if tokens(txt) <= TOK_LIM:
        return [txt]
    chunks, buf, buf_tok = [], [], 0
    for para in txt.splitlines(keepends=True):
        n = tokens(para)
        if buf_tok + n > TOK_LIM:
            chunks.append("".join(buf).strip())
            buf, buf_tok = [], 0
        buf.append(para);  buf_tok += n
    if buf: chunks.append("".join(buf).strip())
    return chunks
# ─────────────────────────────────────────────────────

# ───────────  robust single-chunk writer  ────────────
RETRYABLE   = (RemoteProtocolError, ReadError, APIConnectionError)
MAX_RETRY   = 5
DELAY0      = 2.0        # s
TIMEOUT_REQ = 120        # s

def tts_to_path(text: str, path: Path) -> None:
    """Create one chunk-MP3 with retries (thread-safe)."""
    delay = DELAY0
    for attempt in range(1, MAX_RETRY + 1):
        try:
            mp3 = OpenAI().audio.speech.create(   # new client per thread
                model=MODEL, voice=VOICE, input=text,
                instructions=INSTRUCTIONS, timeout=TIMEOUT_REQ,
            )
            path.write_bytes(mp3.content)
            return
        except RETRYABLE as e:
            if attempt == MAX_RETRY:
                raise
            print(f"  ⚠️  {type(e).__name__} (try {attempt}) – {delay}s back-off")
            time.sleep(delay);  delay *= 2
# ─────────────────────────────────────────────────────

def concat_mp3s(files: Sequence[Path], target: Path) -> None:
    audio = AudioSegment.empty()
    for f in files:
        audio += AudioSegment.from_mp3(f)
    audio.export(target, format="mp3")

# ─────────────  per-section workflow  ───────────────
MAX_WORKERS = 8      # ← parallelism knob

def convert_section(src: Path, dst: Path) -> None:
    if dst.exists():
        print(f"⏩  {dst.name} already complete")
        return

    tmp_dir = dst.with_suffix(".chunks")
    tmp_dir.mkdir(exist_ok=True)
    blocks  = chunk_text(src.read_text())

    # schedule TTS for missing chunks
    with cf.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = []
        for idx, block in enumerate(blocks, 1):
            cpath = tmp_dir / f"chunk_{idx:03d}.mp3"
            if cpath.exists():
                continue
            futures.append(pool.submit(tts_to_path, block, cpath))

        # collect any exceptions early
        for fut in cf.as_completed(futures):
            fut.result()   # will re-raise if tts_to_path failed

    # stitch if every chunk exists
    chunk_files = sorted(tmp_dir.glob("chunk_*.mp3"))
    if len(chunk_files) != len(blocks):
        print(f"✗ {src.name}: still missing chunks – resume next run")
        return

    part = dst.with_suffix(".part.mp3")
    concat_mp3s(chunk_files, part)
    part.rename(dst)
    shutil.rmtree(tmp_dir)
    print(f"✓ {src.name} → {dst}")

# ─────────────  batch driver  ───────────────
def run_batch() -> None:
    SRC_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for txt in sorted(SRC_DIR.glob("*.txt")):
        convert_section(txt, OUT_DIR / f"{txt.stem}.mp3")

# ── execute immediately ──
run_batch()
