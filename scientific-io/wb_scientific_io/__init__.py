"""Generic CSV-only distribution; inherited scientific-format readers excluded.

Gunwoo Park directed the new implementation. MIT, as stated in LICENSE.
"""
__version__ = "0.1.1"


class FormatError(ValueError):
    """An explicitly requested file/column contract could not be satisfied."""
