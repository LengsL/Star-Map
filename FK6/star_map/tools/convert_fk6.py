"""Select SI-mode FK6 astrometry, merge Part I and III, and write binary data.

The resulting file deliberately retains only:
RA(ICRS), DE(ICRS), pmRA* (SI mode), pmDE, and Vmag.
No DRA_* / DDE_* position-difference fields or other prediction modes enter
the output.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from catalog.star import FK6Star, HEADER, MAGIC, RECORD  # noqa: E402


PARTS = ("fk6_part1.ecsv", "fk6_part3.ecsv")


def value_or_nan(value: object) -> float:
    """Convert an Astropy cell to float, representing a masked value as NaN."""
    if bool(getattr(value, "mask", False)) or str(value).strip() in {"", "--"}:
        return math.nan
    return float(value)


def read_part(source: Path) -> list[FK6Star]:
    """Read only the FK6 SI-mode fields requested for the binary catalogue."""
    table = Table.read(source, format="ascii.ecsv")
    required = {"RA_ICRS_", "DE_ICRS_", "pmRA_", "pmDE", "Vmag"}
    missing = required - set(table.colnames)
    if missing:
        raise ValueError(f"{source}: missing required columns {sorted(missing)}")

    stars: list[FK6Star] = []
    for row in table:
        # RA_ICRS_/DE_ICRS_ are the ICRS positions. pmRA_ is already pmRA*
        # in SI mode; do not multiply it by cos(Dec) a second time.
        ra_rad = Angle(str(row["RA_ICRS_"]), unit=u.hourangle).radian
        dec_rad = Angle(str(row["DE_ICRS_"]), unit=u.deg).radian
        stars.append(FK6Star(
            ra_rad=ra_rad,
            dec_rad=dec_rad,
            pmra_star_mas_per_year=value_or_nan(row["pmRA_"]),
            pmdec_mas_per_year=value_or_nan(row["pmDE"]),
            vmag=value_or_nan(row["Vmag"]),
        ))
    return stars


def write_binary(stars: list[FK6Star], destination: Path) -> None:
    """Write a portable little-endian binary: header plus 40-byte records."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with temporary.open("wb") as file:
            file.write(HEADER.pack(MAGIC, len(stars), RECORD.size))
            for star in stars:
                file.write(RECORD.pack(star.ra_rad, star.dec_rad,
                                       star.pmra_star_mas_per_year,
                                       star.pmdec_mas_per_year, star.vmag))
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    raw_dir = CODE_ROOT / "data" / "raw"
    stars = [star for filename in PARTS for star in read_part(raw_dir / filename)]
    expected_count = 878 + 3272
    if len(stars) != expected_count:
        raise RuntimeError(f"expected {expected_count} merged records, found {len(stars)}")
    destination = CODE_ROOT / "data" / "fk6_stars.bin"
    write_binary(stars, destination)
    print(f"Merged {len(stars)} FK6 SI-mode records -> {destination}")
    print(f"Record layout: <5d = RA(rad), Dec(rad), pmRA*(mas/yr), pmDE(mas/yr), Vmag ({RECORD.size} bytes)")


if __name__ == "__main__":
    main()
