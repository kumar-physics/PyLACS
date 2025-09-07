# tests/conftest.py
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]     # repo root
SRC  = ROOT / "src"

# Prefer src/ layout if present
if SRC.exists():
    sys.path.insert(0, str(SRC))
else:
    # fallback: if code is directly under repo root (e.g., repo/pylacs/__init__.py)
    sys.path.insert(0, str(ROOT))
