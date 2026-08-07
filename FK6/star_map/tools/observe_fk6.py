"""Batch-propagate the merged FK6 catalogue and convert it to Alt/Az."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
from pathlib import Path
import sys

import numpy as np

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from astronomy.astronomy import julian_date, radec_to_altaz_batch  # noqa: E402
from astronomy.star_motion import propagate_j2000  # noqa: E402
from catalog.fk6_reader import load_arrays  # noqa: E402


def parse_time(value: str | None) -> datetime:
    if value is None:
        return datetime.now(timezone.utc)
    moment = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return moment.replace(tzinfo=timezone.utc) if moment.tzinfo is None else moment.astimezone(timezone.utc)


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch FK6 J2000 proper motion and Alt/Az conversion")
    parser.add_argument("--time", help="ISO 8601 observation time; default: current UTC")
    parser.add_argument("--longitude", type=float, default=121.55, help="east positive, degrees")
    parser.add_argument("--latitude", type=float, default=29.87, help="north positive, degrees")
    args = parser.parse_args()
    if not -90.0 <= args.latitude <= 90.0:
        raise SystemExit("latitude must be between -90 and 90 degrees")

    observation_time = parse_time(args.time)
    jd_utc = julian_date(observation_time)
    elapsed_years = (jd_utc - 2451545.0) / 365.25
    stars = load_arrays(CODE_ROOT / "data" / "fk6_stars.bin")
    ra, dec = propagate_j2000(stars["ra_rad"], stars["dec_rad"],
                              stars["pmra_star_mas_per_year"], stars["pmdec_mas_per_year"],
                              elapsed_years)
    altitude, azimuth = radec_to_altaz_batch(ra, dec, jd_utc,
                                               math.radians(args.latitude), math.radians(args.longitude))
    above_horizon = altitude >= 0.0
    print(f"Processed {len(stars)} FK6 stars at once")
    print(f"Epoch: {observation_time:%Y-%m-%d %H:%M:%S UTC}; Δt = {elapsed_years:.6f} Julian years")
    print(f"Above horizon: {np.count_nonzero(above_horizon)}")
    for index in range(3):
        print(f"#{index + 1}: Alt={math.degrees(altitude[index]):.5f}°, Az={math.degrees(azimuth[index]):.5f}°")


if __name__ == "__main__":
    main()
