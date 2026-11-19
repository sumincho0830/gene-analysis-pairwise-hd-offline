from __future__ import annotations

import hashlib
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


root = Path(__file__).resolve().parent
failures = []
for line in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    expected, relative = line.split("  ", 1)
    target = root / relative
    if not target.is_file():
        failures.append(f"MISSING {relative}")
    elif sha256(target) != expected:
        failures.append(f"HASH {relative}")

if failures:
    print("FAILED")
    print("\n".join(failures))
    sys.exit(1)
print("OK")
