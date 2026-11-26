"""Header-aware CSV inspection and explicitly selected numeric columns.

Generic new adapter. No scientific meaning, units or transformations inferred.
"""
from __future__ import annotations

import csv
from pathlib import Path
from collections.abc import Sequence

import numpy as np

from . import FormatError

__version__ = "0.1.0"


def _records(path, delimiter, encoding):
    """Read a rectangular CSV and retain physical line numbers for each nonblank record."""
    if not isinstance(delimiter, str) or len(delimiter) != 1:
        raise FormatError("CSV delimiter must be exactly one character")
    path = Path(path)
    rows, blank_records = [], 0
    try:
        with path.open(encoding=encoding, newline="") as handle:
            reader = csv.reader(handle, delimiter=delimiter, strict=True)
            header = next(reader, None)
            if header is None:
                raise FormatError(f"{path}: CSV is empty; a header is required")
            header = [name.strip() for name in header]
            if not header or any(not name for name in header) or len(set(header)) != len(header):
                raise FormatError(f"{path}: CSV header names must be nonempty and unique after stripping")
            for row in reader:
                if not row:
                    blank_records += 1
                    continue
                if len(row) != len(header):
                    raise FormatError(f"{path}:{reader.line_num}: expected {len(header)} fields; got {len(row)}")
                rows.append((reader.line_num, row))
    except csv.Error as error:
        raise FormatError(f"{path}: malformed CSV: {error}") from error
    return path, header, rows, blank_records


def _metadata(path, header, rows, blanks, delimiter, encoding):
    """Count empty, text and non-finite cells per column and report the actual shape."""
    summaries = {}
    for index, name in enumerate(header):
        missing = invalid = nonfinite = 0
        for _, row in rows:
            value = row[index].strip()
            if not value:
                missing += 1
                continue
            try:
                number = float(value)
            except ValueError:
                invalid += 1
            else:
                nonfinite += int(not np.isfinite(number))
        summaries[name] = {"numeric": bool(rows) and not (missing or invalid or nonfinite),
                           "missing_count": missing, "invalid_numeric_count": invalid,
                           "nonfinite_count": nonfinite}
    return {"path": str(path), "format": "csv", "delimiter": delimiter, "encoding": encoding,
            "header": header, "row_count": len(rows), "column_count": len(header),
            "shape": [len(rows), len(header)], "blank_record_count": blanks,
            "columns": summaries, "unit_inference": False}


def inspect_csv(path, *, delimiter=",", encoding="utf-8") -> dict:
    """Inspect actual header/rectangular rows; report missing/text/nonfinite cells.

    Unrelated textual columns are reported, not rejected. Duplicate/empty
    stripped headers and malformed row widths fail before mapping columns.
    """
    path, header, rows, blanks = _records(path, delimiter, encoding)
    return _metadata(path, header, rows, blanks, delimiter, encoding)


def read_csv(path, *, columns: Sequence[str], delimiter=",", encoding="utf-8", missing="error") -> dict:
    """Return selected finite numeric columns plus inspection metadata.

    Only missing='error' is supported: no row dropping, imputation or implicit
    transformation. Requested columns preserve their order and equal lengths.
    Raw unrelated text is neither converted nor treated as a numeric failure.
    """
    if missing != "error":
        raise FormatError("CSV supports missing='error' only; no implicit data transformation")
    if isinstance(columns, (str, bytes)) or not isinstance(columns, Sequence) or not columns:
        raise FormatError("columns must be a nonempty sequence of distinct header names")
    if any(not isinstance(name, str) or not name for name in columns) or len(set(columns)) != len(columns):
        raise FormatError("columns must contain distinct nonempty strings")
    path, header, rows, blanks = _records(path, delimiter, encoding)
    absent = [name for name in columns if name not in header]
    if absent:
        raise FormatError(f"{path}: requested CSV column(s) absent: {', '.join(absent)}")
    if not rows:
        raise FormatError(f"{path}: CSV needs at least one data row")
    selected = {}
    for name in columns:
        index = header.index(name)
        values = []
        for number, row in rows:
            value = row[index].strip()
            if not value:
                raise FormatError(f"{path}:{number}: missing selected value in {name!r}")
            try:
                numeric = float(value)
            except ValueError as error:
                raise FormatError(f"{path}:{number}: invalid numeric value in selected column {name!r}") from error
            if not np.isfinite(numeric):
                raise FormatError(f"{path}:{number}: nonfinite value in selected column {name!r}")
            values.append(numeric)
        selected[name] = np.asarray(values, dtype=np.float64)
    metadata = _metadata(path, header, rows, blanks, delimiter, encoding)
    metadata.update(selected_columns=list(columns), selected_shape=[len(rows), len(columns)], missing_policy=missing)
    return {"columns": selected, "metadata": metadata}
