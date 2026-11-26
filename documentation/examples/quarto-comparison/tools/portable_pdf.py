"""Rewrite Chromium's file URLs after printing a local Quarto slide deck.

Based on the portable PDF link handling used in the wb-client documentation.
"""

from pathlib import Path, PurePosixPath
import argparse
from html.parser import HTMLParser
from urllib.parse import urljoin, urlsplit

from pypdf import PdfWriter
from pypdf.generic import ArrayObject, NameObject, TextStringObject


class DeckLinks(HTMLParser):
    """Collect declared document links and author metadata from the HTML deck."""

    def __init__(self):
        super().__init__()
        self.links = set()
        self.author = None

    def handle_starttag(self, tag, attributes):
        fields = dict(attributes)
        if tag == "meta" and fields.get("name") == "author" and fields.get("content"):
            self.author = fields["content"]
        href = fields.get("href")
        if tag != "a" or not href or href.startswith("#"):
            return
        uri = urlsplit(href)
        if uri.scheme == "https" and uri.hostname and not (uri.username or uri.password):
            self.links.add(href)
            return
        path = PurePosixPath(uri.path)
        if uri.scheme or uri.netloc or path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Expected a portable local deck link: {href!r}")
        self.links.add(href)


def restore_links(html, original, destination):
    """Preserve declared targets and remove Quarto's browser-only line links."""
    deck = DeckLinks()
    deck.feed(html.read_text(encoding="utf-8"))
    declared = {urljoin(html.as_uri(), href): href for href in deck.links}
    runtime = urlsplit(html.as_uri())
    writer = PdfWriter(clone_from=original)
    if deck.author:
        writer.add_metadata({"/Author": deck.author})
    observed = set()
    for page in writer.pages:
        retained = ArrayObject()
        for reference in page.get("/Annots", []):
            annotation = reference.get_object()
            action = annotation.get("/A")
            if annotation.get("/Subtype") == "/Link" and action and action.get("/S") == "/URI":
                uri = str(action["/URI"])
                href = declared.get(uri)
                if href is None:
                    target = urlsplit(uri)
                    if target.scheme == runtime.scheme and target.path == runtime.path and target.query == "print-pdf":
                        continue
                    raise ValueError(f"Unexpected PDF link: {uri!r}")
                action[NameObject("/URI")] = TextStringObject(href)
                observed.add(href)
            retained.append(reference)
        if "/Annots" in page:
            page[NameObject("/Annots")] = retained
    if observed != deck.links:
        raise ValueError(f"Missing PDF links: {sorted(deck.links - observed)}")
    with destination.open("xb") as output:
        writer.write(output)


def main():
    """Select the exact HTML, unmodified browser PDF and fresh portable output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path)
    parser.add_argument("browser_pdf", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    restore_links(args.html.resolve(strict=True), args.browser_pdf, args.destination)


if __name__ == "__main__":
    main()
