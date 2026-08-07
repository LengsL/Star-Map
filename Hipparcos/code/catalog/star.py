"""Build and read the compact Hipparcos star catalogue used by the map.

Data source: CDS/VizieR catalogue I/239, ``hip_main.dat``.  The source
file is a fixed-width / pipe-delimited copy of the ESA Hipparcos Main
Catalogue.  This module uses its RAdeg, DEdeg and Hpmag fields, whose frame
is ICRS at epoch J1991.25.

Run directly to download (when needed), build ``stars.bin``, then print the
first ten binary records::

    python catalog/star.py
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import gzip
import math
from pathlib import Path
import ssl
import struct
import urllib.request
from typing import Iterator

import certifi


HIPPARCOS_URL = "https://cdsarc.cds.unistra.fr/ftp/cats/I/239/hip_main.dat"
MAGIC = b"STARSv1\0"  # 8 bytes
HEADER = struct.Struct("<8sII")  # magic, star count, record size
RECORD = struct.Struct("<Ifff")  # HIP, RA radians, Dec radians, Hpmag


@dataclass(frozen=True)
class Star:
    """One display-ready Hipparcos record."""

    hip: int
    ra_rad: float
    dec_rad: float
    magnitude: float

    @property
    def ra_deg(self) -> float:
        return math.degrees(self.ra_rad)

    @property
    def dec_deg(self) -> float:
        return math.degrees(self.dec_rad)


def download_hipparcos(destination: Path, *, force: bool = False) -> Path:
    """Download the official CDS Hipparcos Main Catalogue."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.stat().st_size > 0 and not force:
        return destination
    request = urllib.request.Request(
        HIPPARCOS_URL,
        headers={"User-Agent": "tianguangsuo-starmap/1.0 (Hipparcos catalogue builder)"},
    )
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        # Use certifi's current CA bundle.  This preserves certificate
        # verification on Windows installations whose OpenSSL store is empty.
        context = ssl.create_default_context(cafile=certifi.where())
        with urllib.request.urlopen(request, timeout=90, context=context) as response, temporary.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def _number(value: str) -> float:
    """Parse VizieR numeric cells and reject blank/unknown values."""
    value = value.strip()
    if not value or value in {"?", "---"}:
        raise ValueError("missing numeric field")
    return float(value)


def _ra_hms_to_degrees(value: str) -> float:
    """Convert Hipparcos H3 right ascension (hh mm ss) to degrees."""
    hours, minutes, seconds = value.split()
    return 15.0 * (float(hours) + float(minutes) / 60.0 + float(seconds) / 3600.0)


def _dec_dms_to_degrees(value: str) -> float:
    """Convert Hipparcos H4 declination (signed dd mm ss) to degrees."""
    degrees, minutes, seconds = value.split()
    sign = -1.0 if degrees.startswith("-") else 1.0
    return sign * (abs(float(degrees)) + float(minutes) / 60.0 + float(seconds) / 3600.0)


def iter_hipparcos(source: Path) -> Iterator[Star]:
    """Yield valid HIP/RAdeg/DEdeg/Hpmag rows from ``hip_main.dat``.

    CDS uses pipes as unambiguous field separators.  Relevant field positions
    are H1 (HIP), H8 (RAdeg), H9 (DEdeg) and H44 (Hpmag).  The V magnitude
    (H5) is intentionally not substituted: records with no Hpmag are omitted
    so the resulting magnitude scale is uniform.
    """
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, mode="rt", encoding="ascii", errors="strict") as catalogue:
        for line_number, row in enumerate(catalogue, start=1):
            fields = row.rstrip("\n").split("|")
            if len(fields) <= 44:
                raise ValueError(f"unexpected I/239 row layout at line {line_number}")
            try:
                hip = int(fields[1])
                hp_mag = _number(fields[44])
            except ValueError:
                continue
            try:
                ra_deg = _number(fields[8])
                dec_deg = _number(fields[9])
            except ValueError:
                # A few catalogue rows lack machine-readable degrees but retain
                # the printed H3/H4 sexagesimal positions.  Preserve them so
                # every Modern IAU constellation endpoint can be connected.
                try:
                    ra_deg = _ra_hms_to_degrees(fields[3])
                    dec_deg = _dec_dms_to_degrees(fields[4])
                except ValueError:
                    continue
            if not (0.0 <= ra_deg < 360.0 and -90.0 <= dec_deg <= 90.0):
                continue
            yield Star(hip, math.radians(ra_deg), math.radians(dec_deg), hp_mag)


def write_stars_bin(stars: Iterator[Star], destination: Path) -> int:
    """Write a versioned little-endian star binary and return its record count."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    count = 0
    try:
        with temporary.open("wb") as output:
            output.write(HEADER.pack(MAGIC, 0, RECORD.size))
            for star in stars:
                output.write(RECORD.pack(star.hip, star.ra_rad, star.dec_rad, star.magnitude))
                count += 1
            output.seek(0)
            output.write(HEADER.pack(MAGIC, count, RECORD.size))
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return count


def read_stars_bin(source: Path) -> Iterator[Star]:
    """Read and validate the binary file created by :func:`write_stars_bin`."""
    with source.open("rb") as catalogue:
        header = catalogue.read(HEADER.size)
        if len(header) != HEADER.size:
            raise ValueError(f"{source}: incomplete stars.bin header")
        magic, count, record_size = HEADER.unpack(header)
        if magic != MAGIC or record_size != RECORD.size:
            raise ValueError(f"{source}: unsupported stars.bin format")
        expected_size = HEADER.size + count * RECORD.size
        if source.stat().st_size != expected_size:
            raise ValueError(f"{source}: size does not match header count")
        for _ in range(count):
            record = catalogue.read(RECORD.size)
            if len(record) != RECORD.size:
                raise ValueError(f"{source}: truncated star record")
            hip, ra_rad, dec_rad, magnitude = RECORD.unpack(record)
            yield Star(hip, ra_rad, dec_rad, magnitude)


def build_catalogue(hipparcos_gz: Path, stars_bin: Path, *, redownload: bool = False) -> int:
    """Ensure the source exists and rebuild the compact binary catalogue."""
    download_hipparcos(hipparcos_gz, force=redownload)
    return write_stars_bin(iter_hipparcos(hipparcos_gz), stars_bin)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1] / "data"
    parser = argparse.ArgumentParser(description="下载并转换 Hipparcos 主星表")
    parser.add_argument("--source", type=Path, default=root / "hip_main.dat")
    parser.add_argument("--output", type=Path, default=root / "stars.bin")
    parser.add_argument("--redownload", action="store_true", help="重新下载原始主星表")
    parser.add_argument("--print-only", action="store_true", help="仅读取既有 stars.bin 并打印前十条")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.print_only:
        if not args.output.is_file():
            raise SystemExit(f"找不到二进制星表：{args.output}")
    else:
        total = build_catalogue(args.source, args.output, redownload=args.redownload)
        print(f"已生成 {args.output}：{total} 颗恒星，每条 {RECORD.size} 字节")
    print("HIP       RA (deg)      Dec (deg)    Hpmag")
    for star in list(read_stars_bin(args.output))[:10]:
        print(f"{star.hip:6d}  {star.ra_deg:12.7f}  {star.dec_deg:12.7f}  {star.magnitude:6.2f}")


if __name__ == "__main__":
    main()
