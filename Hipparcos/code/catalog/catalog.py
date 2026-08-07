"""Compact bright-star catalogue for the visual layer (J2000, radians)."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class Star:
    name: str
    ra_hours: float
    dec_deg: float
    magnitude: float
    constellation: str

    @property
    def ra(self) -> float:
        return math.radians(self.ra_hours * 15.0)

    @property
    def dec(self) -> float:
        return math.radians(self.dec_deg)


STARS = (
    Star("Sirius", 6.7525, -16.7161, -1.46, "Canis Major"),
    Star("Canopus", 6.3992, -52.6957, -0.74, "Carina"),
    Star("Arcturus", 14.2610, 19.1824, -0.05, "Bootes"),
    Star("Vega", 18.6156, 38.7837, 0.03, "Lyra"),
    Star("Capella", 5.2782, 45.9980, 0.08, "Auriga"),
    Star("Rigel", 5.2423, -8.2016, 0.13, "Orion"),
    Star("Procyon", 7.6550, 5.2250, 0.38, "Canis Minor"),
    Star("Betelgeuse", 5.9195, 7.4071, 0.42, "Orion"),
    Star("Achernar", 1.6286, -57.2368, 0.46, "Eridanus"),
    Star("Hadar", 14.0637, -60.3730, 0.61, "Centaurus"),
    Star("Altair", 19.8464, 8.8683, 0.77, "Aquila"),
    Star("Aldebaran", 4.5987, 16.5093, 0.85, "Taurus"),
    Star("Antares", 16.4901, -26.4320, 0.96, "Scorpius"),
    Star("Spica", 13.4199, -11.1613, 0.98, "Virgo"),
    Star("Pollux", 7.7553, 28.0262, 1.14, "Gemini"),
    Star("Deneb", 20.6905, 45.2803, 1.25, "Cygnus"),
    Star("Regulus", 10.1395, 11.9672, 1.35, "Leo"),
    Star("Fomalhaut", 22.9608, -29.6222, 1.16, "Piscis Austrinus"),
    Star("Dubhe", 11.0621, 61.7508, 1.79, "Ursa Major"),
    Star("Merak", 11.0307, 56.3824, 2.37, "Ursa Major"),
    Star("Phecda", 11.8972, 53.6948, 2.44, "Ursa Major"),
    Star("Megrez", 12.2570, 57.0326, 3.31, "Ursa Major"),
    Star("Alioth", 12.9005, 55.9598, 1.76, "Ursa Major"),
    Star("Mizar", 13.3987, 54.9254, 2.23, "Ursa Major"),
    Star("Alkaid", 13.7923, 49.3133, 1.85, "Ursa Major"),
    Star("Alnitak", 5.6793, -1.9426, 1.74, "Orion"),
    Star("Alnilam", 5.6036, -1.2019, 1.69, "Orion"),
    Star("Mintaka", 5.5334, -0.2991, 2.23, "Orion"),
    Star("Saiph", 5.7959, -9.6696, 2.09, "Orion"),
)

CONSTELLATIONS = (
    ("Orion", ("Betelgeuse", "Alnitak", "Alnilam", "Mintaka", "Rigel", "Saiph", "Betelgeuse")),
    ("Ursa Major", ("Dubhe", "Merak", "Phecda", "Megrez", "Alioth", "Mizar", "Alkaid")),
)
