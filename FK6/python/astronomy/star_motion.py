"""Vectorized FK6 proper-motion propagation from J2000 to an observing epoch."""

from __future__ import annotations

import numpy as np


MAS_TO_RADIAN = np.pi / (180.0 * 3_600_000.0)


def propagate_j2000(ra_rad: np.ndarray, dec_rad: np.ndarray,
                    pmra_star_mas_per_year: np.ndarray,
                    pmdec_mas_per_year: np.ndarray,
                    elapsed_years: float) -> tuple[np.ndarray, np.ndarray]:
    """Propagate all stars using the stored FK6 SI-mode proper motions.

    ``pmRA*`` is already ``mu_alpha cos(delta)``.  It is used directly as
    the coefficient of the right-ascension tangent basis vector, avoiding any
    additional multiplication by ``cos(dec)``.
    """
    ra = np.asarray(ra_rad, dtype=np.float64)
    dec = np.asarray(dec_rad, dtype=np.float64)
    mu_alpha_star = np.nan_to_num(np.asarray(pmra_star_mas_per_year, dtype=np.float64)) * MAS_TO_RADIAN
    mu_delta = np.nan_to_num(np.asarray(pmdec_mas_per_year, dtype=np.float64)) * MAS_TO_RADIAN

    sin_ra, cos_ra = np.sin(ra), np.cos(ra)
    sin_dec, cos_dec = np.sin(dec), np.cos(dec)
    position = np.stack((cos_dec * cos_ra, cos_dec * sin_ra, sin_dec), axis=-1)
    # Tangent bases: e_alpha and e_delta.  pmRA* enters e_alpha unchanged.
    e_alpha = np.stack((-sin_ra, cos_ra, np.zeros_like(ra)), axis=-1)
    e_delta = np.stack((-cos_ra * sin_dec, -sin_ra * sin_dec, cos_dec), axis=-1)
    propagated = position + elapsed_years * (mu_alpha_star[:, None] * e_alpha + mu_delta[:, None] * e_delta)
    propagated /= np.linalg.norm(propagated, axis=1)[:, None]

    propagated_ra = np.arctan2(propagated[:, 1], propagated[:, 0]) % (2.0 * np.pi)
    propagated_dec = np.arctan2(propagated[:, 2], np.hypot(propagated[:, 0], propagated[:, 1]))
    return propagated_ra, propagated_dec


def precess_j2000(ra_rad, dec_rad, jd):
    t=(float(jd)-2451545.0)/36525.0
    zeta=np.radians((2306.2181*t+0.30188*t*t+0.017998*t**3)/3600.0); z=np.radians((2306.2181*t+1.09468*t*t+0.018203*t**3)/3600.0); theta=np.radians((2004.3109*t-0.42665*t*t-0.041833*t**3)/3600.0)
    sd,cd=np.sin(dec_rad),np.cos(dec_rad); x,y,q=cd*np.cos(ra_rad),cd*np.sin(ra_rad),sd
    a,b,c=np.cos(theta)*x-np.sin(theta)*q,y,np.sin(theta)*x+np.cos(theta)*q
    u,v=np.cos(zeta)*a-np.sin(zeta)*b,np.sin(zeta)*a+np.cos(zeta)*b
    x2,y2=np.cos(z)*u-np.sin(z)*v,np.sin(z)*u+np.cos(z)*v
    return np.arctan2(y2,x2)%(2*np.pi),np.arctan2(c,np.hypot(x2,y2))
