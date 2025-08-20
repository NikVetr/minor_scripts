#!/usr/bin/env python3
# combine_sections_auto.py  – autodetect output name & chapter order

from pathlib import Path
import re
from pydub import AudioSegment
from mutagen.id3 import ID3, CHAP, CTOC, TIT2, ID3NoHeaderError

# ─────────────── edit just these two if needed ────────────────
VOICE = "shimmer"
MODEL = "gpt-4o-mini-tts"

#SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_creatures_great_and_small/sections")
SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/nothing_much_happens/sections/")
# ──────────────────────────────────────────────────────────────

# input & output locations
IN_DIR   = SECTIONS_DIR / VOICE / MODEL / "mp3"
PROJECT  = SECTIONS_DIR.parent                      # eg …/all_creatures_great_and_small
OUT_MP3  = PROJECT / f"{PROJECT.name}_complete.mp3" # eg all_creatures_great_and_small_complete.mp3

# ---------- gather numbered files ------------------------------------------------
num_re = re.compile(r"^(\d+)-.+\.mp3$", re.I)

def numbered_files(path: Path):
    for f in path.glob("*.mp3"):
        m = num_re.match(f.name)
        if m:
            yield int(m.group(1)), f

files_sorted = [f for _, f in sorted(numbered_files(IN_DIR))]
if not files_sorted:
    raise SystemExit(f"No numbered mp3s found in {IN_DIR}")

# ---------- concatenate while recording chapter offsets -------------------------
combined  = AudioSegment.empty()
starts_ms, labels = [], []

def pretty(stem: str) -> str:
    text = stem.split("-", 1)[1]                # drop leading number-
    return text.replace("_", " ").replace("-", " ").title()

for f in files_sorted:
    starts_ms.append(len(combined))
    combined += AudioSegment.from_mp3(f)
    labels.append(pretty(f.stem))

combined.export(OUT_MP3, format="mp3")
print(f"✓ stitched {len(files_sorted)} files → {OUT_MP3}")

# ---------- embed ID3 v2.4 chapter markers -------------------------------------
try:
    tags = ID3(OUT_MP3)
except ID3NoHeaderError:
    tags = ID3()

chap_ids = []
for ix, (title, start) in enumerate(zip(labels, starts_ms), 1):
    end = starts_ms[ix] if ix < len(starts_ms) else len(combined)
    cid = f"ch{ix:02d}"
    chap_ids.append(cid)
    tags.add(
        CHAP(
            element_id=cid,
            start_time=start,
            end_time=end,
            start_offset=0,
            end_offset=0,
            sub_frames=[TIT2(encoding=3, text=[title])],
        )
    )

tags.add(
    CTOC(
        element_id="toc",
        flags=0,
        child_element_ids=chap_ids,
        sub_frames=[TIT2(encoding=3, text=["Table of Contents"])],
    )
)
tags.save(OUT_MP3, v2_version=4)
print("✓ chapters embedded")
