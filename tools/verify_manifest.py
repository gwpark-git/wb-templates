#!/usr/bin/env python3
"""Check this candidate's bytes/membership without the originating workspace.

Integrity checking does not authenticate a replaced manifest or grant rights.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import stat

MANIFEST = "provenance/DISTRIBUTION_MANIFEST.json"


def verify(root):
    root = Path(root).absolute()
    for path in (*root.parents, root):
        if path.is_symlink():
            raise ValueError("candidate cannot use a symlink route")
    actual = set()
    for current, dirs, files in os.walk(root, followlinks=False):
        for name in dirs + files:
            path = Path(current, name)
            mode = path.lstat().st_mode
            if stat.S_ISLNK(mode) or not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                raise ValueError(f"link/special file: {path.relative_to(root)}")
            if name == ".git":
                if Path(current) != root:
                    raise ValueError("nested Git metadata")
                continue
            if stat.S_ISREG(mode):
                actual.add(path.relative_to(root).as_posix())
        if Path(current) == root:
            dirs[:] = [name for name in dirs if name != ".git"]
    if MANIFEST not in actual:
        raise ValueError("manifest missing")
    manifest = json.loads((root / MANIFEST).read_text(encoding="utf-8"))
    schema = manifest.get("schema")
    if schema not in ("wb.templates.public-distribution.v1", "wb.templates.public-distribution.v2"):
        raise ValueError("unknown manifest schema")
    if schema.endswith(".v2") and set(manifest) != {"schema", "entries", "file_count", "payload_bytes"}:
        raise ValueError("unexpected manifest fields")
    expected, total = {MANIFEST}, 0
    for item in manifest["entries"]:
        if schema.endswith(".v2") and set(item) != {"path", "bytes", "sha256"}:
            raise ValueError("unexpected payload fields")
        name = item["path"]
        path = PurePosixPath(name)
        if (not isinstance(name, str) or str(path) != name or path.is_absolute()
                or ".." in path.parts or ".git" in path.parts or name in expected
                or name not in actual):
            raise ValueError("invalid, duplicate or missing payload path")
        expected.add(name)
        data = (root / name).read_bytes()
        if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise ValueError(f"payload differs: {name}")
        total += len(data)
    if actual != expected:
        raise ValueError(f"membership differs: {sorted(actual ^ expected)}")
    if len(expected) - 1 != manifest["file_count"] or total != manifest["payload_bytes"]:
        raise ValueError("manifest totals differ")
    return {"integrity": "VERIFIED", "file_count": len(expected)-1,
            "payload_bytes": total, "publication_approval": False}


if __name__ == "__main__":
    try:
        print(json.dumps(verify(Path(__file__).absolute().parents[1]), indent=2))
    except (ValueError, KeyError, TypeError, OSError) as error:
        raise SystemExit(f"Verification failed: {error}")
