"""Build a portable student source archive from this checkout; no network calls.

Run: python scripts/package_students.py
The archive includes current working files, not only the last commit.
"""
from pathlib import Path
import argparse
import json
import subprocess
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]


def student_files():
    roots = ["heff", "tests", "examples", "notebooks", "scripts", "docs", "results"]
    allowed = {".py", ".toml", ".md", ".ipynb", ".json", ".csv", ".npz", ".png", ".pdf"}
    files = [ROOT / name for name in ("README.md", "AGENTS.md", "pyproject.toml", ".gitignore", "LICENSE")
             if (ROOT / name).is_file()]
    for name in roots:
        for path in (ROOT / name).rglob("*"):
            relative = path.relative_to(ROOT)
            if not path.is_file() or path.is_symlink() or path.suffix not in allowed:
                continue
            if any(part in {"__pycache__", ".ipynb_checkpoints", ".cache"} for part in relative.parts):
                continue
            # Full papers, book excerpts, and thesis extracts are not teaching dependencies.
            if relative.parts[:2] in (("docs", "thesis-text"), ("docs", "briefs")):
                continue
            if relative.parts[:2] == ("docs", "lit") and path.suffix != ".md":
                continue
            files.append(path)
    return sorted(files)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "dist/heff-students.zip")
    args = parser.parse_args()
    files = student_files()
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        commit, dirty = None, None
    manifest = {"source_commit": commit, "working_tree_dirty": dirty,
                "contents": "Files from the working tree at packaging time",
                "omitted": ["Git metadata", "environment/build/cache files", "coordination checkpoints",
                            "development briefs", "raw literature and thesis extracts"],
                "files": [p.relative_to(ROOT).as_posix() for p in files]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(args.output, "w", ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, "heff/" + path.relative_to(ROOT).as_posix())
        archive.writestr("heff/STUDENT-PACKAGE.json", json.dumps(manifest, indent=2))
    print(f"Created {args.output.resolve()} ({len(files)} files, {args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
