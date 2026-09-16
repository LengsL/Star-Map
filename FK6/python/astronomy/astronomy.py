"""Vectorized RA/Dec to local Alt/Az conversion."""

from __future__ import annotations

from datetime import datetime, timezone
import math

import numpy as np


def julian_date(moment: datetime) -> float:
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc).timestamp() / 86400.0 + 2440587.5


def local_sidereal_time(jd_utc: float, longitude_rad: float) -> float:
    days = jd_utc - 2451545.0
    centuries = days / 36525.0
    gmst_deg = (280.46061837 + 360.98564736629 * days
                + 0.000387933 * centuries * centuries
                - centuries * centuries * centuries / 38710000.0)
    return math.radians(gmst_deg % 360.0) + longitude_rad


def radec_to_altaz_batch(ra_rad: np.ndarray, dec_rad: np.ndarray, jd_utc: float,
                         latitude_rad: float, longitude_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """Convert complete RA/Dec arrays to altitude and north-through-east azimuth."""
    ra = np.asarray(ra_rad, dtype=np.float64)
    dec = np.asarray(dec_rad, dtype=np.float64)
    hour_angle = (local_sidereal_time(jd_utc, longitude_rad) - ra + np.pi) % (2.0 * np.pi) - np.pi
    sin_altitude = (np.sin(dec) * math.sin(latitude_rad)
                    + np.cos(dec) * math.cos(latitude_rad) * np.cos(hour_angle))
    altitude = np.arcsin(np.clip(sin_altitude, -1.0, 1.0))
    cos_altitude = np.maximum(np.cos(altitude), 1e-15)
    sin_azimuth = -np.sin(hour_angle) * np.cos(dec) / cos_altitude
    cos_azimuth = ((np.sin(dec) - np.sin(altitude) * math.sin(latitude_rad))
                   / (cos_altitude * math.cos(latitude_rad)))
    azimuth = np.arctan2(sin_azimuth, cos_azimuth) % (2.0 * np.pi)
    return altitude, azimuth
