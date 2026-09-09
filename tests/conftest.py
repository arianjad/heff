"""Put the repo root on sys.path so tests can import the `scripts/` helpers.

`heff` itself is installed, so tests of the package need nothing here. The
`scripts/` directory is not a package and is not installed; without this,
`python -m pytest` works (it adds the CWD) but a bare `pytest` fails at
collection with "No module named 'scripts'".
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
