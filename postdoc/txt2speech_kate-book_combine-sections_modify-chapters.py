#!/usr/bin/env python3
"""combine_sections_auto_v2.py — chapter trim, gain‑match, **true dual‑silence fade**, stitch, tag

Bug‑fix 2025‑07‑07‑c
────────────────────
*Fade‑out still missing?*  This version reinstates an explicit
`apply_fades()` that calls `AudioSegment.fade()` twice with a user‑set
`FADE_DB` (depth) and *guarantees* the second fade is applied on top of
all prior processing.

Pipeline
────────
1. Trim `CHOP_MS` from the start of each chapter.
2. Gain‑match every chapter to the loudness of chapter 1.
3. `apply_fades()` — linear fade from `FADE_DB` → 0 dB at start, and
   from 0 dB → `FADE_DB` at end.  Default `FADE_DB = -80 dB`, which is
   effectively silence but still measurable.
4. Export processed chapters, print progress.
5. Stitch them, embed ID3 chapters.
"""

from pathlib import Path
import re
from pydub import AudioSegment
from mutagen.id3 import ID3, CHAP, CTOC, TIT2, ID3NoHeaderError

# ─────────────────── CONFIGURATION ───────────────────────────
VOICE  = "shimmer"
MODEL  = "gpt-4o-mini-tts"
#SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/nothing_much_happens/sections/")
#SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_creatures_great_and_small/sections/")
#SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_things_bright_and_beautiful/sections/")
SECTIONS_DIR = Path("/Users/nikgvetr/Documents/Documents - nikolai/all_things_wise_and_wonderful/sections/")

CHOP_MS   = 0     # trim this many ms from the start of every chapter
FADE_MS   = 4000     # fade duration (ms) for both in & out
FADE_DB   = -20.0    # dBFS at fade extremum (≤ –80 ≈ silence)
EXPORT_BITRATE = "128k"
# ─────────────────────────────────────────────────────────────

# ───────────────── PATHS ─────────────────────────────────────
IN_DIR   = SECTIONS_DIR / VOICE / MODEL / "mp3"
PROJECT  = SECTIONS_DIR.parent
PROC_DIR = PROJECT / "processed_chapters"
OUT_MP3  = PROJECT / f"{PROJECT.name}_complete.mp3"
PROC_DIR.mkdir(exist_ok=True)

num_re = re.compile(r"^(\d+)-.+\.mp3$", re.I)

# ───────────────── HELPERS ───────────────────────────────────

def numbered_files(path: Path):
    for f in path.glob("*.mp3"):
        m = num_re.match(f.name)
        if m:
            yield int(m.group(1)), f


def pretty(stem: str) -> str:
    text = stem.split("-", 1)[1]
    return text.replace("_", " ").replace("-", " ").title()


def match_gain(seg: AudioSegment, target_dBFS: float) -> AudioSegment:
    return seg.apply_gain(target_dBFS - seg.dBFS)


def apply_fades(seg: AudioSegment) -> AudioSegment:
    """Fade from `FADE_DB` → 0 dB at start and 0 dB → `FADE_DB` at end."""
    if FADE_MS <= 0 or len(seg) < FADE_MS * 2:
        return seg  # skip if clip too short or fade disabled

    # fade‑in
    seg = seg.fade(start=0,                duration=FADE_MS,
                    from_gain=FADE_DB,     to_gain=0.0)
    # fade‑out
    seg = seg.fade(start=len(seg) - FADE_MS, duration=FADE_MS,
                    from_gain=0.0,         to_gain=FADE_DB)
    return seg

# ───────────────── MAIN ──────────────────────────────────────

files_sorted = [p for _, p in sorted(numbered_files(IN_DIR))]
if not files_sorted:
    raise SystemExit(f"No numbered mp3s found in {IN_DIR}")
print(f"Found {len(files_sorted)} chapter files in {IN_DIR}")

# 1️⃣  target loudness from first chapter (post‑trim)
first_seg = AudioSegment.from_mp3(files_sorted[0])
first_seg = first_seg[CHOP_MS:] if len(first_seg) > CHOP_MS else first_seg
TARGET_dBFS = first_seg.dBFS
print(f"Target dBFS (chapter 1): {TARGET_dBFS:.1f} dBFS\n")

# 2️⃣  per‑chapter processing pass
processed_paths = []
for ix, src in enumerate(files_sorted, 1):
    seg = AudioSegment.from_mp3(src)

    # trim chapter heading
    seg = seg[CHOP_MS:] if len(seg) > CHOP_MS else seg

    # gain‑match to reference loudness
    seg = match_gain(seg, TARGET_dBFS)

    # apply dual fade
    seg = apply_fades(seg)

    out_path = PROC_DIR / src.name
    seg.export(out_path, format="mp3", bitrate=EXPORT_BITRATE)
    processed_paths.append(out_path)
    print(f"[{ix:>3}/{len(files_sorted)}] processed → {out_path.name}  ({len(seg)//1000}s)")

print("\n✓ individual chapters processed\n")

# 3️⃣  concatenate processed chapters
combined   = AudioSegment.empty()
starts_ms, labels = [], []
for ix, p in enumerate(processed_paths, 1):
    starts_ms.append(len(combined))
    combined += AudioSegment.from_mp3(p)
    labels.append(pretty(p.stem))
    print(f"[{ix:>3}/{len(processed_paths)}] stitched {p.name}")

combined.export(OUT_MP3, format="mp3")
print(f"\n✓ wrote combined file → {OUT_MP3}\n")

# 4️⃣  embed ID3 v2.4 chapter markers
try:
    tags = ID3(OUT_MP3)
except ID3NoHeaderError:
    tags = ID3()

chap_ids = []
for ix, (title, start) in enumerate(zip(labels, starts_ms), 1):
    end = starts_ms[ix] if ix < len(starts_ms) else len(combined)
    cid = f"ch{ix:02d}"
    chap_ids.append(cid)
    tags.add(CHAP(element_id=cid, start_time=start, end_time=end, start_offset=0,
                  end_offset=0, sub_frames=[TIT2(encoding=3, text=[title])]))

tags.add(CTOC(element_id="toc", flags=0, child_element_ids=chap_ids,
              sub_frames=[TIT2(encoding=3, text=["Table of Contents"]) ] ))

tags.save(OUT_MP3, v2_version=4)
print("✓ chapters embedded (ID3 v2.4)")
