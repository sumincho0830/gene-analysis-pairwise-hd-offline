from __future__ import annotations

import hashlib
from pathlib import Path
import sys


root = Path(__file__).resolve().parent
failures = []
for line in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    expected, relative = line.split("  ", 1)
    target = root / relative
    if not target.is_file():
        failures.append(f"MISSING {relative}")
        continue
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    if digest != expected:
        failures.append(f"HASH {relative}")

if failures:
    print("FAILED")
    print("\n".join(failures))
    sys.exit(1)
print("OK")
