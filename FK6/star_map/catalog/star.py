"""Binary record definitions for the merged FK6 catalogue."""

from __future__ import annotations

from dataclasses import dataclass
import struct


# Header: magic, record count, bytes per record.
MAGIC = b"FK6BIN1\0"
HEADER = struct.Struct("<8sII")
# Record: RA (radian), Dec (radian), pmRA* (mas/year), pmDec (mas/year), Vmag.
RECORD = struct.Struct("<5d")


@dataclass(frozen=True)
class FK6Star:
    ra_rad: float
    dec_rad: float
    pmra_star_mas_per_year: float
    pmdec_mas_per_year: float
    vmag: float
