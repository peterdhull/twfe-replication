"""Run figure replication from this checkout without installing its package."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from didfigures.__main__ import main

if __name__ == "__main__":
    main()
