#!/usr/bin/env python3
"""Build the complete static wb-templates reader from the public source distribution."""

import argparse
from hashlib import sha256
from html import escape
from html.parser import HTMLParser
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import re
import stat
import subprocess
import tempfile
from urllib.parse import quote, unquote, urlsplit, urlunsplit

DOCUMENTATION_KIT = "documentation/examples/quarto-comparison"
PYTHON_GUIDE_EXAMPLE = "documentation/examples/python-package"
PYTHON_GUIDE_TEMPLATE = "documentation/templates/python-package-guide"
PYTHON_GUIDE_PREVIEW = f"{PYTHON_GUIDE_EXAMPLE}/preview.html"
COLLECTION_GUIDE = "documentation/collection-guide"
COLLECTION_GUIDE_FILES = tuple(f"{COLLECTION_GUIDE}/{name}" for name in (
    "_quarto.yml", "report.qmd", "_guide.md", "_reference.md",
))
COLLECTION_GUIDE_REPORT = f"{COLLECTION_GUIDE}/report.html"
GENERATED_READERS = frozenset({COLLECTION_GUIDE_REPORT, PYTHON_GUIDE_PREVIEW})
READER_RENDER_SOURCES = frozenset({
    *COLLECTION_GUIDE_FILES,
    f"{DOCUMENTATION_KIT}/_publication-info.md",
    f"{PYTHON_GUIDE_EXAMPLE}/guide.md",
    f"{PYTHON_GUIDE_TEMPLATE}/reader.html",
    f"{PYTHON_GUIDE_TEMPLATE}/render.sh",
})
BASIC_READER = "documentation/templates/basic-reader.html"
MANIFEST = "provenance/DISTRIBUTION_MANIFEST.json"
DOC_PROVENANCE = "provenance/DOCUMENTATION_PROVENANCE.json"
RENDER_PROVENANCE = f"{DOCUMENTATION_KIT}/RENDER_PROVENANCE.json"


PARTIALS = {
    BASIC_READER,
    f"{DOCUMENTATION_KIT}/title-slide.html",
    f"{DOCUMENTATION_KIT}/slide-navigation.html",
    f"{PYTHON_GUIDE_TEMPLATE}/reader.html",
}
HREF = re.compile(r'''(?<![\w-])href\s*=\s*(?:"[^"]*"|'[^']*'|[^\s>]+)''', re.I)


def url_relative(target, page):
    return quote(posixpath.relpath(target, posixpath.dirname(page) or "."))


def render(quarto, template, source, data, page, routes):
    """Feed one selected source to a nonexecuting Pandoc reader."""
    org = source.endswith(".org")
    source_format = "org" if org else "markdown+gfm_auto_identifiers-raw_html-raw_attribute"
    if org and re.search(rb"(?im)^\s*#\+(?:include:|setupfile:|begin_export\b)", data):
        raise ValueError(f"{source}: Org includes and export blocks are outside this static build.")
    language = {"README.de.md": "de", "README.ko.md": "ko"}.get(source, "en")
    heading = re.search(rb"(?im)^(?:#\+title:|title:)\s*([^\n]+)$", data)
    if heading is None:
        heading = re.search(rb"(?m)^[#*]{1,6} +([^\n]+)$", data)
    label = heading[1].decode().split(" {#", 1)[0] if heading else PurePosixPath(source).stem
    metadata = {
        "lang": language, "pagetitle": label,
        "home-url": url_relative("index.html", page),
        "documentation-url": url_relative(routes["documentation/README.md"], page),
        "source-url": url_relative(source, page), "source-label": "Org" if org else "Markdown",
    }
    command = [str(quarto), "pandoc", "--sandbox", f"--from={source_format}", "--to=html5",
               "--standalone", "--section-divs", "--toc", "--toc-depth=2", "--no-highlight",
               "--mathml", f"--template={template}"]
    command.extend(f"--metadata={key}:{value}" for key, value in metadata.items())
    result = subprocess.run(command, input=data, capture_output=True, timeout=30)
    if result.returncode or result.stderr:
        raise ValueError(f"{source}: {result.stderr.decode().strip()}")
    return result.stdout.decode()


class Page(HTMLParser):
    """Record only parsed HTML routes, leaving scripts and embedded data untouched."""

    def __init__(self, text):
        super().__init__()
        self.text = text
        self.ids = set()
        self.links = []
        self.line_offsets = [0, *(match.end() for match in re.finditer("\n", text))]
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        line, column = self.getpos()
        self.ids.update(attrs[k] for k in ("id", "name") if attrs.get(k))
        for field in ("href", "src", "poster", "data"):
            if field in attrs and attrs[field] is not None and (field != "data" or tag == "object"):
                self.links.append((tag, field, attrs, line, self.line_offsets[line - 1] + column,
                                   self.get_starttag_text()))


def resolve(page, url):
    """Resolve a local URL against its page inside a portable distribution."""
    parsed = urlsplit(url)
    if parsed.scheme.lower() in {"http", "https", "mailto", "data"} or (parsed.netloc and not parsed.scheme):
        return None
    raw = unquote(parsed.path)
    if parsed.scheme or raw.startswith("/") or "\\" in raw:
        raise ValueError(f"Absolute or unsupported local route: {url}")
    raw, _, locator = raw.partition("::")
    target = posixpath.normpath(posixpath.join(posixpath.dirname(page), raw)) if raw else page
    if target == ".." or target.startswith("../"):
        raise ValueError(f"Route escapes the distribution: {url}")
    fragment = parsed.fragment or (locator[1:] if locator.startswith("#") else "")
    return parsed, target, fragment, locator


def revise_page(name, page, pages, routes, files):
    """Route document links to HTML while retaining explicit edit/source links."""
    edits = []
    for tag, field, attrs, line, offset, raw in page.links:
        url = attrs[field]
        try:
            resolved = resolve(name, url)
            if resolved is None:
                continue
            parts, target, fragment, locator = resolved
            if target not in files:
                raise ValueError(f"Target is not selected: {target}")
            source_link = "data-source-link" in attrs or "download" in attrs
            destination = routes.get(target, target) if tag == "a" and not source_link else target
            if destination in pages and fragment and not fragment.startswith("/"):
                if unquote(fragment) not in pages[destination].ids:
                    raise ValueError(f"Unknown anchor {fragment!r} in {destination}")
            if tag == "a" and field == "href" and (destination != target or locator or "::" in parts.path):
                updated = urlunsplit(("", "", url_relative(destination, name), parts.query, fragment))
                attribute = f'href="{escape(updated, quote=True)}"'
                if locator and not locator.startswith("#"):
                    attribute += f' data-source-location="{escape(locator, quote=True)}"'
                    if "title" not in attrs:
                        attribute += f' title="Source location: {escape(locator, quote=True)}"'
                replacement, count = HREF.subn(lambda _: attribute, raw, count=1)
                if count != 1:
                    raise ValueError(f"Could not rewrite the parsed link: {url}")
                edits.append((offset, offset + len(raw), replacement))
        except ValueError as error:
            raise ValueError(f"{name}:{line}: {error}") from error
    text = page.text
    for start, stop, replacement in reversed(edits):
        text = text[:start] + replacement + text[stop:]
    return text.encode()


def provenance(selected, outputs):
    """Pin the delivered HTMLs while preserving original source/render identities."""
    if DOC_PROVENANCE in selected:
        data = json.loads(selected[DOC_PROVENANCE])
        rows = {row["path"]: row for row in data["files"]
                if row["path"] in selected or row["path"] in outputs}
        for name, payload in selected.items():
            if name in rows or PurePosixPath(name).suffix in {".md", ".org", ".html", ".pdf", ".qmd", ".cff"}:
                rows[name] = {**rows.get(name, {"path": name}), "sha256": sha256(payload).hexdigest()}
        for name, value in outputs.items():
            rows[name] = {**rows.get(name, {"path": name}), "sha256": sha256(value).hexdigest()}
        data["files"] = [rows[name] for name in sorted(rows)]
        yield DOC_PROVENANCE, (json.dumps(data, ensure_ascii=False, indent=2) + "\n").encode()
    if RENDER_PROVENANCE in selected:
        data = json.loads(selected[RENDER_PROVENANCE])
        for name in data["source_sha256"]:
            path = f"{DOCUMENTATION_KIT}/{name}"
            data["source_sha256"][name] = sha256(selected[path]).hexdigest()
        for name in data["output_sha256"]:
            path = f"{DOCUMENTATION_KIT}/{name}"
            if path in outputs:
                data["output_sha256"][name] = sha256(outputs[path]).hexdigest()
        yield RENDER_PROVENANCE, (json.dumps(data, ensure_ascii=False, indent=2) + "\n").encode()


def render_reader_previews(quarto, selected):
    """Build the handbook and Python preview from the exact selected sources."""
    with tempfile.TemporaryDirectory(prefix="wb-reader-previews-") as directory:
        work = Path(directory)
        for name in READER_RENDER_SOURCES:
            source = work / name
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_bytes(selected[name])

        handbook = work / COLLECTION_GUIDE
        subprocess.run([str(quarto), "render", "--no-execute"], cwd=handbook,
                       check=True, timeout=60)
        python_preview = work / PYTHON_GUIDE_PREVIEW
        python_renderer = work / PYTHON_GUIDE_TEMPLATE / "render.sh"
        subprocess.run(["sh", str(python_renderer), "guide.md", python_preview.name],
                       cwd=python_preview.parent, env={**os.environ, "QUARTO": str(quarto)},
                       check=True, timeout=30)
        return {
            COLLECTION_GUIDE_REPORT: Page((handbook / "site/report.html").read_text(encoding="utf-8")),
            PYTHON_GUIDE_PREVIEW: Page(python_preview.read_text(encoding="utf-8")),
        }


def documentation_routes(selected):
    """Assign same-directory HTML destinations to Markdown and Org sources."""
    routes = {name: str(PurePosixPath(name).with_suffix(".html"))
              for name in selected if PurePosixPath(name).suffix in {".md", ".org"}}
    occupied = {name.casefold() for name in selected}
    for route in (*routes.values(), "index.html"):
        if route.casefold() in occupied:
            raise ValueError(f"Documentation HTML route collides with a selected file: {route}")
        occupied.add(route.casefold())
    return routes


def render_site(quarto, template, selected):
    """Render all source documents and resolve links against the complete site."""
    routes = documentation_routes(selected)
    pages = {name: Page(payload.decode()) for name, payload in selected.items()
             if name.endswith(".html") and name not in PARTIALS}
    pages.update(render_reader_previews(quarto, selected))
    for source, page in routes.items():
        pages[page] = Page(render(quarto, template, source, selected[source], page, routes))
    pages["index.html"] = pages[routes["README.md"]]

    files = selected.keys() | pages.keys() | {MANIFEST}
    outputs = {}
    for name, page in pages.items():
        content = revise_page(name, page, pages, routes, files)
        if name not in selected or content != selected[name] or name in GENERATED_READERS:
            outputs[name] = content
    outputs.update(provenance(selected, outputs))
    print(f"Rendered {len(routes)} documents and validated {len(pages)} reader pages.")
    return outputs


def regular_file(path):
    """Read only regular files along paths without symlink components."""
    for component in (*reversed(path.parents), path):
        mode = component.lstat().st_mode
        if stat.S_ISLNK(mode):
            raise ValueError(f"Symlink route refused: {component}")
    if not stat.S_ISREG(mode):
        raise ValueError(f"Expected a regular file: {path}")
    return path.read_bytes()


def read_sources(root):
    """Use the distribution inventory, including local edits to its listed files."""
    manifest = json.loads(regular_file(root / MANIFEST))
    if manifest["schema"] != "wb.templates.public-distribution.v2":
        raise ValueError("A wb-templates source distribution manifest is required.")
    sources = {}
    for entry in manifest["entries"]:
        name = entry["path"]
        path = PurePosixPath(name)
        if (path.is_absolute() or path.as_posix() != name or path == PurePosixPath(".")
                or ".." in path.parts or ".git" in path.parts or "external_references" in path.parts or "\\" in name
                or name == MANIFEST or name in sources):
            raise ValueError(f"Invalid or duplicate inventory path: {name}")
        sources[name] = regular_file(root / name)
    return sources


def build_html(root, quarto, output):
    """Build a self-contained static site from this public source checkout."""
    if output.exists() or output.is_symlink() or output.resolve() != output:
        raise ValueError("Choose a new output directory with no symlink components.")
    if output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError("Choose an output directory outside the checkout.")

    selected = read_sources(root)
    result = {**selected, **render_site(quarto, root / BASIC_READER, selected)}
    manifest = {
        "schema": "wb.templates.public-distribution.v2",
        "entries": [{"path": path, "bytes": len(data), "sha256": sha256(data).hexdigest()}
                    for path, data in sorted(result.items())],
        "file_count": len(result),
        "payload_bytes": sum(map(len, result.values())),
    }
    result[MANIFEST] = (json.dumps(manifest, indent=2, ensure_ascii=False) + "\n").encode()
    output.mkdir()
    for name, data in sorted(result.items()):
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(data)
        source = root / name
        if name in selected:
            path.chmod(stat.S_IMODE(source.stat().st_mode) & 0o777)
    print(f"Open {output / 'index.html'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quarto", required=True, type=Path, help="Path to the Quarto executable")
    parser.add_argument("--output", required=True, type=Path, help="A new directory outside the checkout")
    args = parser.parse_args()
    try:
        output = args.output.expanduser()
        output = output.parent.resolve() / output.name
        build_html(Path(__file__).resolve().parents[1], args.quarto.absolute(), output)
    except (KeyError, OSError, TypeError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(2, f"HTML build refused: {error}\n")


if __name__ == "__main__":
    main()
