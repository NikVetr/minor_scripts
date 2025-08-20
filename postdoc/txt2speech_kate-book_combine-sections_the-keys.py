#!/usr/bin/env python3
# combine_sections.py  – run after your batch script finishes

from pathlib import Path
from pydub import AudioSegment
from mutagen.id3 import ID3, CHAP, CTOC, TIT2, ID3NoHeaderError

# ── locate the finished per-section files ───────────────────────────
VOICE  = "shimmer"
MODEL  = "gpt-4o-mini-tts"

#BASE   = Path("/Users/nikgvetr/Documents/Documents - nikolai/the_keys/sections")
BASE   = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_creatures_great_and_small/sections")
IN_DIR = BASE / VOICE / MODEL / "mp3"
OUT_MP3 = BASE / "the-keys_complete.mp3"

SECTION_ORDER = [
    "the-keys_prologue",
    "the-keys_part-1_gathering",
    "the-keys_part-2_capture-and-escape",
    "the-keys_part-3_when-vardelien-dies",
    "the-keys_part-4_creatures-of-the-mines",
    "the-keys_part-5_battle-of-darien…-yet-again",
    "the-keys_part-6_dissembled-feelings",
    "the-keys_epilogue",
]

# ── 1. concatenate while capturing start times ───────────────────────
combined   = AudioSegment.empty()
start_ms   = []                       # chapter start offsets in milliseconds

for name in SECTION_ORDER:
    mp3_path = IN_DIR / f"{name}.mp3"
    if not mp3_path.exists():
        raise FileNotFoundError(mp3_path)
    start_ms.append(len(combined))
    combined += AudioSegment.from_mp3(mp3_path)

combined.export(OUT_MP3, format="mp3")
print(f"✓ stitched → {OUT_MP3}")

# ── 2. add ID3 chapter frames (v2.4) ─────────────────────────────────
try:
    tags = ID3(OUT_MP3)
except ID3NoHeaderError:
    tags = ID3()                      # create new tag set

chap_ids = []
for ix, (name, start) in enumerate(zip(SECTION_ORDER, start_ms), 1):
    end = start_ms[ix] if ix < len(start_ms) else len(combined)
    cid = f"ch{ix:02d}"
    chap_ids.append(cid)

    tags.add(
        CHAP(
            element_id=cid,
            start_time=start,
            end_time=end,
            start_offset=0,
            end_offset=0,
            sub_frames=[TIT2(encoding=3, text=[name.replace('-', ' ').title()])],
        )
    )

# single table of contents referencing all chapters
tags.add(
    CTOC(
        element_id="toc",
        flags=0,
        child_element_ids=chap_ids,
        sub_frames=[TIT2(encoding=3, text=["Table of Contents"])],
    )
)

tags.save(OUT_MP3, v2_version=4)      # ID3 v2.4 needed for CHAP/CTOC
print("✓ chapters embedded (ID3 v2.4)")
