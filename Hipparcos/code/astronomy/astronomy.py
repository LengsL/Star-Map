"""Convert ICRS right ascension/declination to observer Alt/Az coordinates."""

from __future__ import annotations

from datetime import datetime, timezone
import math

AU_KM = 149_597_870.7
EARTH_EQUATORIAL_RADIUS_KM = 6_378.137


def julian_date(moment: datetime) -> float:
    """Convert a datetime to UTC Julian Date."""
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc).timestamp() / 86400.0 + 2440587.5


def local_sidereal_time(jd_utc: float, longitude_rad: float) -> float:
    """Return local mean sidereal time in radians; east longitude is positive."""
    days = jd_utc - 2451545.0
    centuries = days / 36525.0
    gmst_deg = (280.46061837 + 360.98564736629 * days
                + 0.000387933 * centuries * centuries
                - centuries * centuries * centuries / 38710000.0)
    return math.radians(gmst_deg % 360.0) + longitude_rad


def radec_to_altaz(ra_rad: float, dec_rad: float, jd_utc: float,
                   latitude_rad: float, longitude_rad: float) -> tuple[float, float]:
    """Return altitude and north-through-east azimuth in radians."""
    lst = local_sidereal_time(jd_utc, longitude_rad)
    hour_angle = (lst - ra_rad + math.pi) % (2.0 * math.pi) - math.pi
    sin_altitude = (math.sin(dec_rad) * math.sin(latitude_rad)
                    + math.cos(dec_rad) * math.cos(latitude_rad) * math.cos(hour_angle))
    altitude = math.asin(max(-1.0, min(1.0, sin_altitude)))
    cos_altitude = max(1e-15, math.cos(altitude))
    sin_azimuth = -math.sin(hour_angle) * math.cos(dec_rad) / cos_altitude
    cos_azimuth = ((math.sin(dec_rad) - math.sin(altitude) * math.sin(latitude_rad))
                   / (cos_altitude * math.cos(latitude_rad)))
    azimuth = math.atan2(sin_azimuth, cos_azimuth) % (2.0 * math.pi)
    return altitude, azimuth


def vector_to_radec(vector: tuple[float, float, float]) -> tuple[float, float]:
    """Turn an equatorial Cartesian direction vector into RA/Dec."""
    x, y, z = vector
    return math.atan2(y, x) % (2.0 * math.pi), math.atan2(z, math.hypot(x, y))


def topocentric_vector(geocentric_vector_au: tuple[float, float, float], jd_utc: float,
                       latitude_rad: float, longitude_rad: float) -> tuple[float, float, float]:
    """Approximate observer-topocentric vector in the mean equatorial frame."""
    lst = local_sidereal_time(jd_utc, longitude_rad)
    radius_au = EARTH_EQUATORIAL_RADIUS_KM / AU_KM
    observer = (
        radius_au * math.cos(latitude_rad) * math.cos(lst),
        radius_au * math.cos(latitude_rad) * math.sin(lst),
        radius_au * math.sin(latitude_rad),
    )
    return tuple(target - site for target, site in zip(geocentric_vector_au, observer))
