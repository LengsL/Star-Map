from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from astronomy.astronomy import radec_to_altaz_batch
from astronomy.star_motion import propagate_j2000
from catalog.fk6_reader import iter_stars, load_arrays


class FK6CoordinateTest(unittest.TestCase):
    def test_icrs_angles_and_motion_fields_are_present(self) -> None:
        first = next(iter_stars(ROOT / "data" / "fk6_stars.bin"))
        self.assertTrue(0.0 <= first.ra_rad < 2.0 * math.pi)
        self.assertTrue(-math.pi / 2.0 <= first.dec_rad <= math.pi / 2.0)
        self.assertTrue(math.isfinite(first.pmra_star_mas_per_year))
        self.assertTrue(math.isfinite(first.pmdec_mas_per_year))

    def test_all_stars_can_be_propagated_and_converted(self) -> None:
        stars = load_arrays(ROOT / "data" / "fk6_stars.bin")
        ra, dec = propagate_j2000(stars["ra_rad"], stars["dec_rad"],
                                  stars["pmra_star_mas_per_year"], stars["pmdec_mas_per_year"], 26.0)
        altitude, azimuth = radec_to_altaz_batch(ra, dec, 2460000.5, math.radians(29.87), math.radians(121.55))
        self.assertEqual(len(altitude), 4150)
        self.assertTrue(((altitude >= -math.pi / 2) & (altitude <= math.pi / 2)).all())
        self.assertTrue(((azimuth >= 0.0) & (azimuth < 2.0 * math.pi)).all())


if __name__ == "__main__":
    unittest.main()
