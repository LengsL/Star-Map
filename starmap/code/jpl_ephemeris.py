"""Reader for JPL direct-access binary ephemerides (DE200/DE406).

The implementation follows the public USNO/NOVAS ``eph_manager.c``
algorithm supplied with this project.  It intentionally does not treat the
files as SPICE BSP kernels: they are a different, older JPL binary format.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import struct
from typing import Final

AU_KM: Final = 149_597_870.7


class EphemerisError(RuntimeError):
    """Base error raised while opening or evaluating an ephemeris."""


class DateOutOfRange(EphemerisError):
    """The requested Julian date is outside the file's coverage."""


@dataclass(frozen=True)
class EphemerisHeader:
    start_jd: float
    end_jd: float
    step_days: float
    de_number: int
    au_km: float
    earth_moon_ratio: float
    pointers: tuple[tuple[int, int, int], ...]
    record_length: int


class JPLEphemeris:
    """Evaluate positions in an original JPL DE direct-access file.

    Body indices are NOVAS indices: Mercury=0 … Pluto=8, geocentric Moon=9,
    Sun=10, solar-system barycentre=11, and Earth-Moon barycentre=12.
    Returned vectors are AU and AU/day in the file's mean-equator frame.
    """

    _RECORD_LENGTHS: Final = {200: 6608, 403: 8144, 404: 5824,
                              405: 8144, 406: 5824, 421: 8144}

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._file = self.path.open("rb")
        try:
            self.header = self._read_header()
        except Exception:
            self._file.close()
            raise
        self._cached_record_number: int | None = None
        self._cached_record: tuple[float, ...] | None = None

    def close(self) -> None:
        if not self._file.closed:
            self._file.close()

    def __enter__(self) -> "JPLEphemeris":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _read_header(self) -> EphemerisHeader:
        # The supplied lnxm files are little-endian, IEEE-754 64-bit doubles.
        raw = self._file.read(2856)
        if len(raw) != 2856:
            raise EphemerisError(f"{self.path.name}: header is incomplete")
        offset = 252 + 2400  # title and constant-name text fields
        start, end, step = struct.unpack_from("<3d", raw, offset)
        offset += 24
        _ncon = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        au_km, em_ratio = struct.unpack_from("<2d", raw, offset)
        offset += 16
        ints = struct.unpack_from("<36i", raw, offset)
        # File order is one (start, coefficient-count, subinterval-count)
        # triple per body.  The supplied NOVAS C reader scatters those triples
        # into its C[3][12] array while reading them.
        pointers = tuple((ints[3 * i], ints[3 * i + 1], ints[3 * i + 2])
                         for i in range(12))
        offset += 36 * 4
        de_number = struct.unpack_from("<i", raw, offset)[0]
        record_length = self._RECORD_LENGTHS.get(de_number)
        if record_length is None:
            raise EphemerisError(
                f"{self.path.name}: unsupported DE{de_number}; "
                "this reader supports DE200/403/404/405/406/421")
        # DE406 begins at JD 625360.5 (year −3000), so do not assume a
        # modern-era Julian-date lower bound here.
        if not (0 < start < end < 4_000_000 and step > 0):
            raise EphemerisError(f"{self.path.name}: invalid header or byte order")
        return EphemerisHeader(start, end, step, de_number, au_km, em_ratio,
                               pointers, record_length)

    def state(self, jd_tdb: float, body: int) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Return a raw barycentric state (Moon is geocentric), in AU/AU-day."""
        if not 0 <= body <= 10:
            raise ValueError("raw state body must be in 0..10")
        h = self.header
        if not h.start_jd <= jd_tdb <= h.end_jd:
            raise DateOutOfRange(
                f"JD {jd_tdb:.5f} is outside DE{h.de_number} "
                f"({h.start_jd:.1f}–{h.end_jd:.1f})")
        # Record 1 is header, record 2 constants, record 3 first coefficient set.
        record = int((jd_tdb - h.start_jd) / h.step_days) + 3
        if jd_tdb == h.end_jd:
            record -= 2
        record_start = h.start_jd + (record - 3) * h.step_days
        fraction = (jd_tdb - record_start) / h.step_days
        coefficients = self._get_record(record)
        start, ncf, na = h.pointers[body]
        if start <= 0 or ncf <= 0 or na <= 0:
            raise EphemerisError(f"DE{h.de_number} has no coefficients for body {body}")
        pos_km, vel_km_per_day = self._interpolate(
            coefficients, start - 1, ncf, na, fraction, h.step_days)
        return (tuple(v / h.au_km for v in pos_km),
                tuple(v / h.au_km for v in vel_km_per_day))

    def position(self, jd_tdb: float, target: int, center: int = 2) -> tuple[float, float, float]:
        """Return target relative to center; default is Earth-centered."""
        if target == center:
            return (0.0, 0.0, 0.0)
        if target == 11:
            target_pos = (0.0, 0.0, 0.0)
        elif target == 12:
            target_pos = self.state(jd_tdb, 2)[0]
        else:
            target_pos = self.state(jd_tdb, target)[0]
        if center == 11:
            center_pos = (0.0, 0.0, 0.0)
        elif center == 12:
            center_pos = self.state(jd_tdb, 2)[0]
        else:
            center_pos = self.state(jd_tdb, center)[0]

        moon = self.state(jd_tdb, 9)[0] if (target == 2 or center == 2) else None
        earth_bary = self.state(jd_tdb, 2)[0] if (target in (9, 12) or center in (9, 12)) else None
        if target == 2:
            target_pos = _sub(target_pos, _scale(moon, 1.0 / (1.0 + self.header.earth_moon_ratio)))
        elif target == 9:
            target_pos = _add(_sub(earth_bary, _scale(target_pos, 1.0 / (1.0 + self.header.earth_moon_ratio))), target_pos)
        if center == 2:
            center_pos = _sub(center_pos, _scale(moon, 1.0 / (1.0 + self.header.earth_moon_ratio)))
        elif center == 9:
            center_pos = _add(_sub(earth_bary, _scale(center_pos, 1.0 / (1.0 + self.header.earth_moon_ratio))), center_pos)
        return _sub(target_pos, center_pos)

    def _get_record(self, number: int) -> tuple[float, ...]:
        if number != self._cached_record_number:
            length = self.header.record_length
            self._file.seek((number - 1) * length)
            raw = self._file.read(length)
            if len(raw) != length:
                raise EphemerisError(f"could not read coefficient record {number}")
            self._cached_record = struct.unpack(f"<{length // 8}d", raw)
            self._cached_record_number = number
        assert self._cached_record is not None
        return self._cached_record

    @staticmethod
    def _interpolate(buf: tuple[float, ...], start: int, ncf: int, na: int,
                     fraction: float, interval_days: float) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        sub = min(na - 1, int(na * fraction))
        tc = 2.0 * ((na * fraction) % 1.0) - 1.0
        polynomials = [1.0, tc]
        derivs = [0.0, 1.0]
        for i in range(2, ncf):
            polynomials.append(2.0 * tc * polynomials[-1] - polynomials[-2])
            derivs.append(2.0 * tc * derivs[-1] + 2.0 * polynomials[-2] - derivs[-2])
        position, velocity = [], []
        factor = 2.0 * na / interval_days
        for component in range(3):
            first = start + sub * 3 * ncf + component * ncf
            coeff = buf[first:first + ncf]
            position.append(sum(p * c for p, c in zip(polynomials, coeff)))
            velocity.append(factor * sum(d * c for d, c in zip(derivs, coeff)))
        return tuple(position), tuple(velocity)


def _add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def _sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def _scale(a: tuple[float, float, float] | None, value: float) -> tuple[float, float, float]:
    assert a is not None
    return tuple(x * value for x in a)  # type: ignore[return-value]
