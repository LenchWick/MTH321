"""Verify the imported October archive without rerunning Monte Carlo experiments."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[2]
manifest = json.loads((Path(__file__).with_name("archive_sha256.json")).read_text())
for name, expected in manifest.items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise AssertionError(f"Archive checksum mismatch: {name}")
print(f"Verified {len(manifest)} archived Project 1 files.")
