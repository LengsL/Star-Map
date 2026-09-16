"""Validated reader for ``data/fk6_stars.bin``."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np

from catalog.star import FK6Star, HEADER, MAGIC, RECORD


ARRAY_DTYPE = np.dtype([
    ("ra_rad", "<f8"),
    ("dec_rad", "<f8"),
    ("pmra_star_mas_per_year", "<f8"),
    ("pmdec_mas_per_year", "<f8"),
    ("vmag", "<f8"),
])


def read_header(source: Path) -> int:
    """Validate the FK6 binary file and return its number of star records."""
    with source.open("rb") as file:
        raw = file.read(HEADER.size)
    if len(raw) != HEADER.size:
        raise ValueError(f"{source}: incomplete FK6 header")
    magic, count, record_size = HEADER.unpack(raw)
    if magic != MAGIC or record_size != RECORD.size:
        raise ValueError(f"{source}: unsupported FK6 binary format")
    expected_size = HEADER.size + count * RECORD.size
    if source.stat().st_size != expected_size:
        raise ValueError(f"{source}: binary size does not match header")
    return count


def iter_stars(source: Path) -> Iterator[FK6Star]:
    """Yield stars in the merged FK6 binary in Part I then Part III order."""
    count = read_header(source)
    with source.open("rb") as file:
        file.seek(HEADER.size)
        for _ in range(count):
            raw = file.read(RECORD.size)
            if len(raw) != RECORD.size:
                raise ValueError(f"{source}: truncated FK6 record")
            yield FK6Star(*RECORD.unpack(raw))


def load_stars(source: Path) -> list[FK6Star]:
    """Load all FK6 records into memory."""
    return list(iter_stars(source))


def load_arrays(source: Path) -> np.ndarray:
    """Read every FK6 record at once into a NumPy structured array."""
    count = read_header(source)
    with source.open("rb") as file:
        file.seek(HEADER.size)
        stars = np.fromfile(file, dtype=ARRAY_DTYPE, count=count)
    if len(stars) != count:
        raise ValueError(f"{source}: truncated FK6 records")
    return stars
