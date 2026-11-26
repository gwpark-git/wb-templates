#!/usr/bin/env python3
"""Copy a template or collection without local runtime/Git state.

This tool performs a local source copy only; it installs and executes nothing.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import stat


EXCLUDED = frozenset({
    ".git", ".work", "external_references", "build", "runs", ".venv", "venv", "__pycache__",
    ".pytest_cache", ".mypy_cache", ".ruff_cache", ".DS_Store",
    "wb-templates-public",
})


def excluded(name: str) -> bool:
    return (name in EXCLUDED or name.endswith((".pyc", ".pyo", ".egg-info"))
            or (name.startswith("#") and name.endswith("#")))


def omitted_name(directory: Path, name: str) -> bool:
    return excluded(name) or (directory.name == "publication" and name == "history")


def inspect_source(source: Path) -> None:
    """Reject links before copying so a source tree cannot import outside data."""
    for directory, directories, files in os.walk(source, followlinks=False):
        directories[:] = [name for name in directories if not omitted_name(Path(directory), name)]
        for name in directories + [name for name in files if not excluded(name)]:
            path = Path(directory, name)
            mode = path.lstat().st_mode
            if stat.S_ISLNK(mode):
                raise ValueError(f"Source symlink is not distributable: {path.relative_to(source)}")
            if not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                raise ValueError(f"Source special file is not distributable: {path.relative_to(source)}")


def export(root: Path, template: str | None, destination: Path) -> Path:
    root = root.resolve()
    catalog = json.loads((root / "templates.json").read_text(encoding="utf-8"))
    selected = {item["id"]: item for item in catalog["templates"] if item.get("kind") != "documentation"}
    if template is not None and template not in selected:
        raise ValueError(f"Unknown template {template!r}; choose {', '.join(selected)}")
    item = selected[template] if template is not None else {}
    source_path = root if template is None else root / item["path"]
    if source_path.is_symlink():
        raise ValueError("Catalog source must not be a symlink")
    source = source_path.resolve()
    if not source.is_relative_to(root):
        raise ValueError("Catalog source escapes the collection")
    if not source.is_dir():
        raise ValueError(f"Template source is missing: {source.name}")
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError("Destination already exists; choose a new directory")
    if destination.resolve().is_relative_to(source):
        raise ValueError("Destination must be outside the selected source tree")
    if not (source / "LICENSE").is_file() and template is not None:
        raise ValueError("Independent template must carry its own LICENSE")
    inspect_source(source)
    # A detached mBLA copy carries the collection MIT grant under its established name.
    license_copy = item.get("shared_license_copy")
    if license_copy is not None:
        if license_copy != "LICENSE.NEW":
            raise ValueError("Unsupported shared license destination")
        grant = root / "LICENSE"
        if grant.is_symlink() or not grant.is_file() or (source / license_copy).exists():
            raise ValueError("The shared license needs one regular source and a fresh destination")
    shutil.copytree(source, destination, ignore=lambda _path, names: [n for n in names if omitted_name(Path(_path), n)])
    if license_copy is not None:
        shutil.copy2(grant, destination / license_copy)
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", help="Template ID from templates.json; omit for the collection")
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    try:
        destination = export(Path(__file__).resolve().parents[1], args.template, args.destination)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as error:
        parser.exit(2, f"Export refused: {error}\n")
    print(f"Source copied to {destination}; dependencies were not installed and no model was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
