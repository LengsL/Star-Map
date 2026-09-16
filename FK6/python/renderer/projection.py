from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime, timezone
import math
from pathlib import Path
import sys
import time

import numpy as np
import pygame
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers

# Never access the network during startup or coordinate refresh.  Astropy will
# use its bundled/local IERS table; this is sufficient for the interactive map
# and avoids downloading a refreshed table on every launch.
iers.conf.auto_download = False
iers.conf.auto_max_age = None

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from astronomy.astronomy import julian_date, local_sidereal_time, radec_to_altaz_batch  # noqa: E402
from astronomy.star_motion import propagate_j2000  # noqa: E402
from catalog.fk6_reader import load_arrays  # noqa: E402


# Visual palette: muted colours preserve dark adaptation and separate layers.
BACKGROUND = (3, 6, 12)
GRID = (38, 48, 62)
GRID_LABEL = (120, 133, 152)
HORIZON = (88, 157, 181)
HOUR_ANGLE = (120, 94, 136)
EQUATOR = (73, 116, 164)
ECLIPTIC = (175, 144, 85)
CONSTELLATION = (105, 85, 125)
NORTH_POLE = (93, 163, 192)
HUD_BACKGROUND = (13, 20, 31)
HUD_BORDER = (76, 123, 151)
HUD_TEXT = (210, 220, 234)
HUD_ALPHA = 205
SIDEBAR_WIDTH = 350
SIDEBAR_GAP = 28
MAX_PAN = 1.25  # Horizon-radius units; applies to the complete view transform.
DECLINATION_GRID_DEGREES = (-60, -30, 0, 30, 60)
HOUR_ANGLE_GRID_HOURS = tuple(range(-12, 12))
MIN_WINDOW_WIDTH, MIN_WINDOW_HEIGHT = 760, 560
TELESCOPE = (75, 220, 140)

# J2000 RA/Dec (degrees) for familiar asterism lines.  Keeping these compact
# coordinates locally makes the FK6 renderer self-contained: it never reads a
MAJOR_ASTERISMS = {
    "Ori": ((88.79, 7.41), (81.28, 6.35), (83.00, -0.30), (84.05, -1.20),
            (85.19, -1.94), (86.94, -9.67), (78.63, -8.20)),
    "UMa": ((165.93, 61.75), (165.46, 56.38), (178.46, 53.69), (183.86, 57.03),
            (193.51, 55.96), (200.98, 54.93), (206.89, 49.31)),
    "Cas": ((2.29, 59.15), (10.13, 56.54), (21.45, 60.24), (28.60, 63.67), (50.00, 56.54)),
    "Cyg": ((310.36, 45.28), (305.56, 40.26), (292.68, 27.96), (305.25, 33.97), (296.24, 45.13)),
    "Sco": ((247.35, -26.43), (252.17, -43.00), (254.66, -37.10), (263.40, -37.10),
            (264.33, -42.00), (265.62, -39.03)),
    "Leo": ((152.09, 11.97), (154.99, 19.84), (154.17, 23.42), (168.53, 20.52),
            (168.56, 15.43), (177.26, 14.57)),
    "Lyr": ((279.23, 38.78), (283.82, 43.95), (281.19, 37.61), (284.74, 32.69)),
    "Aql": ((297.70, 8.87), (291.37, 3.11), (286.56, -4.88), (305.56, 12.44)),
    "Sgr": ((276.04, -34.38), (283.82, -26.30), (286.74, -27.67), (290.97, -29.88),
            (283.54, -30.42)),
    "Peg": ((345.94, 28.08), (346.19, 15.21), (326.05, 9.88), (340.75, 30.22)),
    "Tau": ((68.98, 16.51), (65.73, 17.54), (64.95, 15.63), (81.57, 28.61)),
}
CONSTELLATION_NAMES = {
    "Ori": "猎户座", "UMa": "大熊座", "Cas": "仙后座", "Cyg": "天鹅座",
    "Sco": "天蝎座", "Leo": "狮子座", "Lyr": "天琴座", "Aql": "天鹰座",
    "Sgr": "人马座", "Peg": "飞马座", "Tau": "金牛座",
}
BODY_STYLES = {
    "Sun": ((237, 198, 119), 7), "Moon": ((192, 204, 219), 6),
    "Mercury": ((156, 150, 142), 3), "Venus": ((219, 184, 140), 4),
    "Mars": ((193, 111, 92), 4), "Jupiter": ((205, 172, 141), 5),
    "Saturn": ((196, 185, 135), 5), "Uranus": ((112, 190, 201), 3),
    "Neptune": ((102, 133, 199), 3), "Pluto": ((151, 137, 126), 2),
}



    
@dataclass(frozen=True)
class SolarSystemBody:
    name: str
    altitude_rad: float
    azimuth_rad: float
    colour: tuple[int, int, int]
    radius: int


@dataclass(frozen=True)
class MoonPhaseInfo:
    illumination: float
    waxing: bool


@dataclass(frozen=True)
class TelescopePointing:
    altitude_deg: float
    azimuth_deg: float


@dataclass(frozen=True)
class StarHit:
    index: int
    x: int
    y: int
    radius: int


@dataclass
class LabelPlacer:
    """Place map labels only when their screen rectangles do not overlap."""

    rectangles: list[pygame.Rect]

    def __init__(self) -> None:
        self.rectangles = []

    def place(self, surface: pygame.Surface, text: pygame.Surface,
              position: tuple[int, int], *, centered: bool = False) -> bool:
        rectangle = text.get_rect(center=position) if centered else text.get_rect(topleft=position)
        if any(rectangle.colliderect(existing) for existing in self.rectangles):
            return False
        surface.blit(text, rectangle)
        self.rectangles.append(rectangle)
        return True

    def place_near(self, surface: pygame.Surface, text: pygame.Surface,
                   point: tuple[int, int], offsets: tuple[tuple[int, int], ...]) -> bool:
        for dx, dy in offsets:
            if self.place(surface, text, (point[0] + dx, point[1] + dy)):
                return True
        return False


def ecliptic_radec(samples: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """J2000 ecliptic expressed as equatorial RA/Dec samples."""
    longitude = np.linspace(0.0, 2.0 * math.pi, samples)
    obliquity = math.radians(23.4392911)
    return (np.arctan2(np.sin(longitude) * math.cos(obliquity), np.cos(longitude)) % (2.0 * math.pi),
            np.arcsin(np.sin(longitude) * math.sin(obliquity)))


def load_major_constellation_lines() -> list[tuple[str, np.ndarray, np.ndarray]]:
    """Return self-contained J2000 major-constellation line coordinates."""
    return [(name, np.radians(np.asarray(points)[:, 0]), np.radians(np.asarray(points)[:, 1]))
            for name, points in MAJOR_ASTERISMS.items()]


def project_altaz(altitude: np.ndarray | float, azimuth: np.ndarray | float,
                  centre: tuple[float, float], radius: float, zoom: float,
                  pan: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """Project cached horizon coordinates in a zenithal-equidistant chart."""
    altitude_array = np.asarray(altitude, dtype=np.float64)
    azimuth_array = np.asarray(azimuth, dtype=np.float64)
    distance = radius * zoom * (math.pi / 2.0 - altitude_array) / (math.pi / 2.0)
    return (centre[0] + pan[0] + distance * np.sin(azimuth_array),
            centre[1] + pan[1] - distance * np.cos(azimuth_array))


@dataclass
class FK6SkyState:
    source: Path
    latitude_deg: float = 29.87
    longitude_deg: float = 121.55
    magnitude_limit: float | None = None
    motion_epoch: datetime | None = None

    telescope_pointing: tuple[float, float] | None = None


    def set_telescope_pointing(
        self,
        altitude_deg:float,
        azimuth_deg:float
    ) -> None:
        if not -90.0 <= altitude_deg <= 90.0:
            raise ValueError("Altitude must be between -90.0 to 90.0 degree")

        if not 0.0 <= azimuth_deg < 360.0:
            raise ValueError("Azimuth must be between 0.0 to 360.0 degree")
        self.telescope_pointing = TelescopePointing(altitude_deg, azimuth_deg)




    
    def __post_init__(self) -> None:

        self.source = Path(self.source)
        self.stars = load_arrays(self.source)  # Read FK6 exactly once.
        self.latitude_rad = math.radians(self.latitude_deg)
        self.longitude_rad = math.radians(self.longitude_deg)
        finite_magnitudes = self.stars["vmag"][np.isfinite(self.stars["vmag"])]
        if not len(finite_magnitudes):
            raise ValueError("FK6 catalogue contains no finite Vmag values")
        self.catalog_mag_min = math.floor(float(finite_magnitudes.min()))
        self.catalog_mag_max = math.ceil(float(finite_magnitudes.max()))
        self.magnitude_limit = (float(self.catalog_mag_max) if self.magnitude_limit is None
                                else self._clamp_magnitude(self.magnitude_limit))
        self.zoom = 1.0
        # In horizon-radius units.  It is constrained to the zoom surplus so
        # the enlarged sky always covers the fixed circular coordinate chart.
        self.pan = (0.0, 0.0)
        self.propagation_count = 0
        self.coordinate_update_count = 0
        self.projection_dirty = True

        self.motion_epoch = self.motion_epoch or datetime.now(timezone.utc)
        if self.motion_epoch.tzinfo is None:
            self.motion_epoch = self.motion_epoch.replace(tzinfo=timezone.utc)
        elapsed_years = (julian_date(self.motion_epoch) - 2451545.0) / 365.25
        # This is the only proper-motion propagation for this program run.
        self.ra_rad, self.dec_rad = propagate_j2000(
            self.stars["ra_rad"], self.stars["dec_rad"],
            self.stars["pmra_star_mas_per_year"], self.stars["pmdec_mas_per_year"],
            elapsed_years,
        )

        self.propagation_count += 1
        self.altitude_rad = np.empty(len(self.stars), dtype=np.float64)
        self.azimuth_rad = np.empty(len(self.stars), dtype=np.float64)
        self.hour_angle_lines: list[tuple[int, np.ndarray, np.ndarray]] = []
        self.equator_radec = (np.linspace(0.0, 2.0 * math.pi, 361), np.zeros(361))
        self.ecliptic_radec = ecliptic_radec()
        self.reference_altaz: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self.north_pole_altaz = (0.0, 0.0)
        self.solar_system: list[SolarSystemBody] = []
        self.layers = {"constellations": True, "equatorial": True,
                       "altaz_grid": True, "hour_angle": True}
        self.constellation_edges = load_major_constellation_lines()
        self.constellation_altaz: list[tuple[str, np.ndarray, np.ndarray]] = []
        self.observation_time = self.motion_epoch
        self.visible_star_hits: list[StarHit] = []
        self.goto_star: int | None = None
        self.selected_star_index: int | None = None
        self.target_goto_rect: pygame.Rect | None = None
        self.moon_phase = MoonPhaseInfo(0.0, True)
        self.update_altaz(self.motion_epoch)

    def _clamp_magnitude(self, value: float) -> float:
        return max(float(self.catalog_mag_min), min(float(self.catalog_mag_max), float(value)))

    def update_altaz(self, moment: datetime) -> None:
        """Refresh cached local coordinates; no J2000 propagation happens here."""
        if moment.tzinfo is None:
            moment = moment.replace(tzinfo=timezone.utc)
        self.observation_time = moment.astimezone(timezone.utc)
        jd = julian_date(self.observation_time)
        self.altitude_rad, self.azimuth_rad = radec_to_altaz_batch(
            self.ra_rad, self.dec_rad, jd, self.latitude_rad, self.longitude_rad)
        self.hour_angle_lines = self._make_hour_angle_lines(jd)
        self.reference_altaz = {
            "equator": radec_to_altaz_batch(*self.equator_radec, jd, self.latitude_rad, self.longitude_rad),
            "ecliptic": radec_to_altaz_batch(*self.ecliptic_radec, jd, self.latitude_rad, self.longitude_rad),
        }
        for declination in DECLINATION_GRID_DEGREES:
            ra = np.linspace(0.0, 2.0 * math.pi, 361)
            dec = np.full_like(ra, math.radians(declination))
            self.reference_altaz[f"dec_{declination:+d}"] = radec_to_altaz_batch(
                ra, dec, jd, self.latitude_rad, self.longitude_rad)
        north_alt, north_az = radec_to_altaz_batch(np.array((0.0,)), np.array((math.pi / 2.0,)),
                                                    jd, self.latitude_rad, self.longitude_rad)
        self.north_pole_altaz = (float(north_alt[0]), float(north_az[0]))
        self.constellation_altaz = []
        for abbreviation, ra, dec in self.constellation_edges:
            altitude, azimuth = radec_to_altaz_batch(ra, dec, jd, self.latitude_rad, self.longitude_rad)
            self.constellation_altaz.append((abbreviation, altitude, azimuth))
        self.solar_system = self._solar_system_positions()
        self.coordinate_update_count += 1
        self.projection_dirty = True

    def _solar_system_positions(self) -> list[SolarSystemBody]:
        """Cache topocentric Sun, Moon and planetary Alt/Az for one second."""
        observing_time = Time(self.observation_time)
        location = EarthLocation.from_geodetic(self.longitude_deg * u.deg, self.latitude_deg * u.deg)
        frame = AltAz(obstime=observing_time, location=location)
        bodies: list[SolarSystemBody] = []
        for name in BODY_STYLES:
            try:
                coordinate = (get_sun(observing_time) if name == "Sun"
                              else get_body(name.lower(), observing_time, location))
                horizon = coordinate.transform_to(frame)
                colour, radius = BODY_STYLES[name]
                bodies.append(SolarSystemBody(name, float(horizon.alt.to_value(u.rad)),
                                              float(horizon.az.to_value(u.rad)), colour, radius))
                if name == "Moon":
                    moon = coordinate.cartesian.xyz.to_value(u.au)
                    sun = get_sun(observing_time).cartesian.xyz.to_value(u.au)
                    elongation = math.acos(float(np.dot(moon, sun) /
                                                 max(np.linalg.norm(moon) * np.linalg.norm(sun), 1e-12)))
                    illumination = (1.0 - math.cos(elongation)) / 2.0
                    # The sign of the ecliptic longitude difference distinguishes waxing/waning.
                    self.moon_phase = MoonPhaseInfo(float(np.clip(illumination, 0.0, 1.0)),
                                                     bool((moon[0] * sun[1] - moon[1] * sun[0]) >= 0.0))
            except Exception:
                # The built-in ephemeris may not offer every distant body in a
                # particular Astropy version; the remaining bodies still draw.
                continue
        return bodies

    def _make_hour_angle_lines(self, jd: float) -> list[tuple[int, np.ndarray, np.ndarray]]:
        """Cache constant-hour-angle lines, sampled by declination."""
        lst = local_sidereal_time(jd, self.longitude_rad)
        dec = np.linspace(-math.pi / 2.0 + 1e-5, math.pi / 2.0 - 1e-5, 181)
        lines = []
        for hour in HOUR_ANGLE_GRID_HOURS:
            ra = np.full_like(dec, (lst - hour * math.pi / 12.0) % (2.0 * math.pi))
            altitude, azimuth = radec_to_altaz_batch(ra, dec, jd, self.latitude_rad, self.longitude_rad)
            lines.append((hour, altitude, azimuth))
        return lines

    # The following controls only mark the projection surface dirty.
    def zoom_by(self, multiplier: float) -> None:
        self.zoom = max(1.0, min(4.0, self.zoom * multiplier))
        self._clamp_pan()
        self.projection_dirty = True

    def _clamp_pan(self) -> None:
        """Keep the complete movable view within a finite navigation range."""
        maximum = MAX_PAN
        distance = math.hypot(*self.pan)
        if distance > maximum and distance > 0.0:
            scale = maximum / distance
            self.pan = (self.pan[0] * scale, self.pan[1] * scale)

    def pan_by(self, dx: float, dy: float) -> None:
        """Move every map layer together in the shared view transform."""
        self.pan = (self.pan[0] + dx, self.pan[1] + dy)
        self._clamp_pan()
        self.projection_dirty = True

    def adjust_magnitude(self, delta: float) -> None:
        self.magnitude_limit = self._clamp_magnitude(self.magnitude_limit + delta)
        self.projection_dirty = True

    def reset_projection(self) -> None:
        self.zoom = 1.0
        self.pan = (0.0, 0.0)
        self.projection_dirty = True

    def goto_selected_target(self) -> bool:
        if self.selected_star_index is None:
            return False

        index = self.selected_star_index
        self.goto_star = index

        self.set_telescope_pointing(math.degrees(float(self.altitude_rad[index])),
                                    math.degrees(float(self.azimuth_rad[index])))
        self.projection_dirty = True

        return True



def _view_pan(state: FK6SkyState, radius: int) -> tuple[float, float]:
    return state.pan[0] * radius, state.pan[1] * radius


def view_geometry(state: FK6SkyState, centre: tuple[int, int], radius: int) -> tuple[tuple[float, float], float]:
    """Screen centre and horizon radius after the shared zoom/pan transform."""
    pan_x, pan_y = _view_pan(state, radius)
    return (centre[0] + pan_x, centre[1] + pan_y), radius * state.zoom


def _inside_view(x: np.ndarray | float, y: np.ndarray | float,
                 view_centre: tuple[float, float], view_radius: float) -> np.ndarray:
    return ((np.asarray(x) - view_centre[0]) ** 2 + (np.asarray(y) - view_centre[1]) ** 2
            <= view_radius ** 2)


def draw_horizon_grid(surface: pygame.Surface, font: pygame.font.Font, state: FK6SkyState,
                      centre: tuple[int, int], radius: int, labels: LabelPlacer) -> None:
    """Draw the horizon grid and labels through the same view transform as stars."""
    view_centre, view_radius = view_geometry(state, centre, radius)
    view_centre_int = (round(view_centre[0]), round(view_centre[1]))
    if state.layers["altaz_grid"]:
        for altitude_deg in (15, 30, 45, 60, 75):
            ring_radius = round(view_radius * (90 - altitude_deg) / 90)
            pygame.draw.circle(surface, GRID, view_centre_int, ring_radius, 1)
            label = font.render(f"Alt {altitude_deg}°", True, GRID_LABEL)
            labels.place(surface, label, (view_centre_int[0] + 5, view_centre_int[1] - ring_radius - 16))
        for azimuth_deg in range(0, 360, 30):
            x, y = project_altaz(0.0, math.radians(azimuth_deg), centre, radius,
                                 state.zoom, _view_pan(state, radius))
            pygame.draw.line(surface, GRID, view_centre_int, (round(float(x)), round(float(y))), 1)
    pygame.draw.circle(surface, HORIZON, view_centre_int, round(view_radius), 2)
    for azimuth_deg, label in ((0, "N"), (90, "E"), (180, "S"), (270, "W")):
        x, y = project_altaz(0.0, math.radians(azimuth_deg), centre, radius,
                             state.zoom, _view_pan(state, radius))
        direction_x, direction_y = float(x) - view_centre[0], float(y) - view_centre[1]
        length = max(1.0, math.hypot(direction_x, direction_y))
        text = font.render(label, True, HORIZON)
        labels.place(surface, text,
                     (round(float(x) + 19 * direction_x / length),
                      round(float(y) + 19 * direction_y / length)), centered=True)
    if state.layers["altaz_grid"]:
        pygame.draw.circle(surface, HORIZON, view_centre_int, 4, 1)
        zenith = font.render("Zenith", True, HORIZON)
        labels.place(surface, zenith, (view_centre_int[0], view_centre_int[1] + 18), centered=True)


def draw_hour_angle_lines(surface: pygame.Surface, state: FK6SkyState,
                          font: pygame.font.Font, centre: tuple[int, int], radius: int,
                          labels: LabelPlacer) -> None:
    """Draw cached hour-angle curves after the horizon coordinate grid."""
    if not state.layers["hour_angle"]:
        return
    view_centre, view_radius = view_geometry(state, centre, radius)
    for hour, altitude, azimuth in state.hour_angle_lines:
        x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom,
                             (state.pan[0] * radius, state.pan[1] * radius))
        inside = (altitude >= 0.0) & _inside_view(x, y, view_centre, view_radius)
        segment: list[tuple[int, int]] = []
        label_position: tuple[int, int] | None = None
        for current_x, current_y, visible in zip(x, y, inside):
            if visible:
                point = (round(float(current_x)), round(float(current_y)))
                segment.append(point)
                label_position = point
            else:
                if len(segment) > 1:
                    pygame.draw.aalines(surface, HOUR_ANGLE, False, segment)
                segment = []
        if len(segment) > 1:
            pygame.draw.aalines(surface, HOUR_ANGLE, False, segment)
        if label_position is not None:
            label = font.render(f"H {hour:+d}h", True, HOUR_ANGLE)
            labels.place(surface, label, label_position, centered=True)


def draw_cached_curve(surface: pygame.Surface, altitude: np.ndarray, azimuth: np.ndarray,
                      colour: tuple[int, int, int], state: FK6SkyState,
                      centre: tuple[int, int], radius: int) -> None:
    """Draw the visible parts of a cached Alt/Az curve inside the horizon."""
    x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom, _view_pan(state, radius))
    view_centre, view_radius = view_geometry(state, centre, radius)
    inside = (altitude >= 0.0) & _inside_view(x, y, view_centre, view_radius)
    segment: list[tuple[int, int]] = []
    for current_x, current_y, visible in zip(x, y, inside):
        if visible:
            segment.append((round(float(current_x)), round(float(current_y))))
        else:
            if len(segment) > 1:
                pygame.draw.aalines(surface, colour, False, segment)
            segment = []
    if len(segment) > 1:
        pygame.draw.aalines(surface, colour, False, segment)


def draw_reference_layers(surface: pygame.Surface, state: FK6SkyState,
                          centre: tuple[int, int], radius: int) -> None:
    """Draw cached equator, ecliptic and locally-defined constellation lines."""
    if state.layers["equatorial"]:
        altitude, azimuth = state.reference_altaz["equator"]
        draw_cached_curve(surface, altitude, azimuth, EQUATOR, state, centre, radius)
        for declination in DECLINATION_GRID_DEGREES:
            altitude, azimuth = state.reference_altaz[f"dec_{declination:+d}"]
            draw_cached_curve(surface, altitude, azimuth, EQUATOR, state, centre, radius)
    if state.layers["equatorial"]:
        altitude, azimuth = state.reference_altaz["ecliptic"]
        draw_cached_curve(surface, altitude, azimuth, ECLIPTIC, state, centre, radius)
    if state.layers["constellations"]:
        for _, altitude, azimuth in state.constellation_altaz:
            draw_cached_curve(surface, altitude, azimuth, CONSTELLATION, state, centre, radius)


def draw_constellation_labels(surface: pygame.Surface, state: FK6SkyState,
                              font: pygame.font.Font, centre: tuple[int, int], radius: int,
                              labels: LabelPlacer) -> None:
    """Place visible constellation names from the same cached Alt/Az lines."""
    if not state.layers["constellations"]:
        return
    view_centre, view_radius = view_geometry(state, centre, radius)
    for abbreviation, altitude, azimuth in state.constellation_altaz:
        x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom, _view_pan(state, radius))
        visible = (altitude >= 0.0) & _inside_view(x, y, view_centre, view_radius)
        if not np.any(visible):
            continue
        # Average the visible projected points, so labels track zoom, pan and
        # the once-per-second Alt/Az refresh without introducing new geometry.
        point = (round(float(np.mean(x[visible]))), round(float(np.mean(y[visible]))))
        label = font.render(CONSTELLATION_NAMES.get(abbreviation, abbreviation), True, CONSTELLATION)
        labels.place_near(surface, label, point,
                          ((7, -17), (7, 7), (-42, -17), (-42, 7), (7, -33), (-42, 23)))


def draw_stars(surface: pygame.Surface, state: FK6SkyState,
               centre: tuple[int, int], radius: int) -> int:
    """Draw FK6 stars; their appearance is determined only by their own Vmag."""
    x, y = project_altaz(state.altitude_rad, state.azimuth_rad, centre, radius, state.zoom,
                         _view_pan(state, radius))
    magnitudes = state.stars["vmag"]
    view_centre, view_radius = view_geometry(state, centre, radius)
    inside = _inside_view(x, y, view_centre, view_radius)
    visible = ((state.altitude_rad >= 0.0) & np.isfinite(magnitudes)
               & (magnitudes <= state.magnitude_limit) & inside)
    state.visible_star_hits = []

    for index in np.flatnonzero(visible):
        brightness, dot_radius = star_appearance(float(magnitudes[index]))
        point = (round(float(x[index])), round(float(y[index])))
        state.visible_star_hits.append(StarHit(int(index), point[0], point[1], max(2, dot_radius)))

        colour = (brightness, min(255, brightness + 4), min(255, brightness + 12))
        if dot_radius == 0:
            surface.set_at(point, colour)
            continue
        if dot_radius >= 3:
            halo_radius = dot_radius + 3
            halo = pygame.Surface((halo_radius * 2 + 2, halo_radius * 2 + 2), pygame.SRCALPHA)
            pygame.draw.circle(halo, (*colour, 30), (halo_radius + 1, halo_radius + 1), halo_radius)
            surface.blit(halo, (point[0] - halo_radius - 1, point[1] - halo_radius - 1))
        pygame.draw.circle(surface, colour, point, dot_radius)
    return int(visible.sum())


def find_hovered_star(state: FK6SkyState, mouse_position: tuple[int, int]) -> StarHit | None:
    """Find the nearest visible star in screen space, independent of sky coordinates."""
    if not state.visible_star_hits:
        return None
    mx, my = mouse_position
    hit = min(state.visible_star_hits, key=lambda item: (item.x - mx) ** 2 + (item.y - my) ** 2)
    limit = max(7, hit.radius + 5)
    return hit if (hit.x - mx) ** 2 + (hit.y - my) ** 2 <= limit * limit else None


def draw_star_tooltip(surface: pygame.Surface, state: FK6SkyState,
                      font: pygame.font.Font, mouse_position: tuple[int, int]) -> None:
    hit = find_hovered_star(state, mouse_position)
    if hit is None:
        return
    i = hit.index
    ra_hours = math.degrees(float(state.ra_rad[i])) / 15.0 % 24.0
    dec_deg = math.degrees(float(state.dec_rad[i]))
    alt_deg = math.degrees(float(state.altitude_rad[i]))
    az_deg = math.degrees(float(state.azimuth_rad[i]))
    lines = (f"FK6 index: {i}", f"Vmag: {state.stars['vmag'][i]:.2f}",
             f"RA: {ra_hours:.3f} h", f"Dec: {dec_deg:+.3f}°",
             f"Alt: {alt_deg:+.2f}°", f"Az: {az_deg:.2f}°")
    rendered = [font.render(line, True, HUD_TEXT) for line in lines]
    width = max(item.get_width() for item in rendered) + 16
    height = sum(item.get_height() for item in rendered) + 12
    x, y = mouse_position[0] + 14, mouse_position[1] + 14
    if x + width > surface.get_width(): x = mouse_position[0] - width - 14
    if y + height > surface.get_height(): y = mouse_position[1] - height - 14
    panel = pygame.Surface((width, height), pygame.SRCALPHA)
    pygame.draw.rect(panel, (*HUD_BACKGROUND, 230), panel.get_rect(), border_radius=5)
    pygame.draw.rect(panel, (*HUD_BORDER, 230), panel.get_rect(), 1, border_radius=5)
    for row, text in enumerate(rendered): panel.blit(text, (8, 6 + row * text.get_height()))
    surface.blit(panel, (x, y))


def _format_ra(hours: float) -> str:
    hours %= 24.0
    hour = int(hours)
    minute_total = (hours - hour) * 60.0
    minute = int(minute_total)
    second = round((minute_total - minute) * 60.0)
    if second == 60:
        second, minute = 0, minute + 1
    return f"{hour:02d}h {minute:02d}m {second:02d}s"


def draw_selected_target_panel(surface: pygame.Surface, state: FK6SkyState,
                               font: pygame.font.Font) -> None:
    index = state.selected_star_index
    if index is None:
        state.target_goto_rect = None
        return
    width, height = 300, 238
    panel = pygame.Rect(18, surface.get_height() - height - 18, width, height)
    card = pygame.Surface(panel.size, pygame.SRCALPHA)
    pygame.draw.rect(card, (*HUD_BACKGROUND, 238), card.get_rect(), border_radius=8)
    pygame.draw.rect(card, (*HUD_BORDER, 235), card.get_rect(), 1, border_radius=8)
    ra_hours = math.degrees(float(state.ra_rad[index])) / 15.0
    dec = math.degrees(float(state.dec_rad[index]))
    alt = math.degrees(float(state.altitude_rad[index]))
    az = math.degrees(float(state.azimuth_rad[index]))
    lines = ("Selected Target", f"index: {index}", f"RA: {_format_ra(ra_hours)}",
             f"Dec: {dec:+07.2f}°", f"Vmag: {float(state.stars['vmag'][index]):.2f}",
             "", f"Current Alt: {alt:5.1f}°", f"Current Az:  {az:5.1f}°")
    for row, line in enumerate(lines):
        card.blit(font.render(line, True, HUD_TEXT if row == 0 else GRID_LABEL), (14, 10 + row * 23))
    button = pygame.Rect(104, height - 42, 92, 29)
    pygame.draw.rect(card, (38, 102, 78), button, border_radius=5)
    pygame.draw.rect(card, TELESCOPE, button, 1, border_radius=5)
    label = font.render("GOTO", True, (235, 255, 240))
    card.blit(label, label.get_rect(center=button.center))
    surface.blit(card, panel.topleft)
    state.target_goto_rect = pygame.Rect(panel.x + button.x, panel.y + button.y,
                                         button.width, button.height)


def draw_object_layers(surface: pygame.Surface, state: FK6SkyState,
                       font: pygame.font.Font, centre: tuple[int, int], radius: int,
                       labels: LabelPlacer) -> None:
    """Draw NCP then Sun, Moon and planets on top of the stellar catalogue."""
    view_centre, view_radius = view_geometry(state, centre, radius)
    if state.layers["equatorial"]:
        altitude, azimuth = state.north_pole_altaz
        if altitude >= 0.0:
            x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom, _view_pan(state, radius))
            point = (round(float(x)), round(float(y)))
            if math.dist(point, view_centre) <= view_radius:
                pygame.draw.circle(surface, NORTH_POLE, point, 6, 1)
                pygame.draw.line(surface, NORTH_POLE, (point[0] - 8, point[1]), (point[0] + 8, point[1]), 1)
                pygame.draw.line(surface, NORTH_POLE, (point[0], point[1] - 8), (point[0], point[1] + 8), 1)
                labels.place_near(surface, font.render("NCP", True, NORTH_POLE), point,
                                  ((9, -17), (9, 7), (-38, -17), (-38, 7)))
    for body in state.solar_system:
        if body.name == "Moon":
            continue
        if body.altitude_rad < 0.0:
            continue
        x, y = project_altaz(body.altitude_rad, body.azimuth_rad, centre, radius, state.zoom,
                             _view_pan(state, radius))
        point = (round(float(x)), round(float(y)))
        if math.dist(point, view_centre) <= view_radius:
            pygame.draw.circle(surface, body.colour, point, body.radius)
            labels.place_near(surface, font.render(body.name, True, body.colour), point,
                              ((body.radius + 4, -10), (body.radius + 4, 7),
                               (-body.radius - 50, -10), (-body.radius - 50, 7)))


def draw_moon_layer(surface: pygame.Surface, state: FK6SkyState,
                    font: pygame.font.Font, centre: tuple[int, int], radius: int,
                    labels: LabelPlacer) -> None:
    """Draw the Moon as a small disc with its illuminated fraction, not a white dot."""
    body = next((item for item in state.solar_system if item.name == "Moon"), None)
    if body is None or body.altitude_rad < 0.0:
        return
    view_centre, view_radius = view_geometry(state, centre, radius)
    x, y = project_altaz(body.altitude_rad, body.azimuth_rad, centre, radius, state.zoom,
                         _view_pan(state, radius))
    point = (round(float(x)), round(float(y)))
    if math.dist(point, view_centre) > view_radius:
        return
    r = max(5, body.radius + 1)
    pygame.draw.circle(surface, (35, 42, 55), point, r)
    pygame.draw.circle(surface, (205, 211, 220), point, r, 1)
    illumination = state.moon_phase.illumination
    if illumination > 0.001:
        # A vertical ellipse gives a readable crescent/gibbous approximation.
        lit_width = max(1, round(2 * r * illumination))
        lit_rect = pygame.Rect(point[0] - r if state.moon_phase.waxing else point[0] + r - lit_width,
                               point[1] - r, lit_width, 2 * r)
        pygame.draw.ellipse(surface, (225, 229, 235), lit_rect)
    labels.place_near(surface, font.render(f"Moon {illumination * 100:.0f}%", True, body.colour), point,
                      ((r + 4, -10), (r + 4, 7), (-r - 65, -10), (-r - 65, 7)))


def star_appearance(magnitude: float) -> tuple[int, int]:
    """Stable Vmag-to-dot mapping, independent of the selected Vmag limit."""
    level = max(0.0, min(1.0, (8.5 - magnitude) / 10.0))
    brightness = round(105 + 150 * (level ** 0.45))
    radius = (4 if magnitude <= 0.0 else 3 if magnitude <= 2.0 else
              2 if magnitude <= 4.0 else 1 if magnitude <= 6.0 else 0)
    return brightness, radius


def render_frame(state: FK6SkyState, width: int, height: int,
                 font: pygame.font.Font, hud_font: pygame.font.Font,
                 mouse_position: tuple[int, int] | None = None) -> tuple[pygame.Surface, int]:
    """Reproject only already-cached Alt/Az coordinates onto a new surface."""
    surface = pygame.Surface((width, height))
    surface.fill(BACKGROUND)
    map_left = min(width - 260, SIDEBAR_WIDTH + SIDEBAR_GAP)
    map_width = width - map_left - 18
    centre = (map_left + map_width // 2, height // 2 + 24)
    radius = max(100, min(map_width // 2 - 20, height // 2 - 62))
    labels = LabelPlacer()
    draw_horizon_grid(surface, font, state, centre, radius, labels)
    draw_hour_angle_lines(surface, state, font, centre, radius, labels)
    draw_reference_layers(surface, state, centre, radius)
    visible = draw_stars(surface, state, centre, radius)
    draw_constellation_labels(surface, state, font, centre, radius, labels)
    draw_object_layers(surface, state, font, centre, radius, labels)
    draw_moon_layer(surface, state, font, centre, radius, labels)
    if state.telescope_pointing is not None:
        pointing = state.telescope_pointing
        altitude = pointing.altitude_deg if hasattr(pointing, "altitude_deg") else pointing[0]
        azimuth = pointing.azimuth_deg if hasattr(pointing, "azimuth_deg") else pointing[1]
        x, y = project_altaz(math.radians(altitude), math.radians(azimuth), centre, radius,
                             state.zoom, _view_pan(state, radius))
        point = (round(float(x)), round(float(y)))
        pygame.draw.circle(surface, TELESCOPE, point, 8, 1)
        pygame.draw.line(surface, TELESCOPE, (point[0] - 11, point[1]), (point[0] + 11, point[1]), 1)
        pygame.draw.line(surface, TELESCOPE, (point[0], point[1] - 11), (point[0], point[1] + 11), 1)

    panel = pygame.Rect(16, 16, SIDEBAR_WIDTH - 32, 158)
    panel_surface = pygame.Surface(panel.size, pygame.SRCALPHA)
    pygame.draw.rect(panel_surface, (*HUD_BACKGROUND, HUD_ALPHA), panel_surface.get_rect(), border_radius=8)
    pygame.draw.rect(panel_surface, (*HUD_BORDER, 225), panel_surface.get_rect(), 1, border_radius=8)
    surface.blit(panel_surface, panel)
    details = (
        "FK6 Horizon Sky Map",
        f"UTC: {state.observation_time:%Y-%m-%d %H:%M:%S}",
        f"Stars: {visible:,} / {len(state.stars):,}   Vmag ≤ {state.magnitude_limit:.1f}",
        "Wheel: zoom   Arrows/WASD: move",
        "1 Constellations   2 Equatorial",
        "3 Alt/Az grid   4 Hour Angle   0 Reset   Esc Quit",
    )
    for line_number, detail in enumerate(details):
        colour = HUD_TEXT if line_number == 0 else GRID_LABEL
        surface.blit(hud_font.render(detail, True, colour), (panel.x + 12, panel.y + 10 + 21 * line_number))
    draw_selected_target_panel(surface, state, font)
    if mouse_position is not None:
        draw_star_tooltip(surface, state, font, mouse_position)
    state.projection_dirty = False
    return surface, visible


def parse_time(value: str | None) -> datetime | None:
    if value is None:
        return None
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return parsed.replace(tzinfo=timezone.utc) if parsed.tzinfo is None else parsed


def apply_event(event: pygame.event.Event, state: FK6SkyState) -> bool:
    """Apply display-only controls. Return ``False`` when the app should quit."""
    if event.type == pygame.MOUSEBUTTONDOWN:
        if event.button == 1:
            if state.target_goto_rect is not None and state.target_goto_rect.collidepoint(event.pos):
                state.goto_selected_target()
            else:
                hit = find_hovered_star(state, event.pos)
                if hit is not None:
                    state.selected_star_index = hit.index
                    state.selected_star = hit.index
                    state.projection_dirty = True

    if event.type == pygame.QUIT:
        return False
    if event.type == pygame.MOUSEWHEEL:
        state.zoom_by(1.15 if event.y > 0 else 1.0 / 1.15)
    if event.type == pygame.KEYDOWN:
        key_actions = {
            pygame.K_UP: (0.0, -0.06), pygame.K_w: (0.0, -0.06),
            pygame.K_DOWN: (0.0, 0.06), pygame.K_s: (0.0, 0.06),
            pygame.K_LEFT: (-0.06, 0.0), pygame.K_a: (-0.06, 0.0),
            pygame.K_RIGHT: (0.06, 0.0), pygame.K_d: (0.06, 0.0),
        }
        if event.key == pygame.K_ESCAPE:
            return False
        if event.key in key_actions:
            state.pan_by(*key_actions[event.key])
        elif event.key == pygame.K_LEFTBRACKET:
            state.adjust_magnitude(-0.25)
        elif event.key == pygame.K_RIGHTBRACKET:
            state.adjust_magnitude(0.25)
        elif event.key == pygame.K_0:
            state.reset_projection()
        elif event.key == pygame.K_1:
            state.layers["constellations"] = not state.layers["constellations"]
            state.projection_dirty = True
        elif event.key == pygame.K_2:
            state.layers["equatorial"] = not state.layers["equatorial"]
            state.projection_dirty = True
        elif event.key == pygame.K_3:
            state.layers["altaz_grid"] = not state.layers["altaz_grid"]
            state.projection_dirty = True
        elif event.key == pygame.K_4:
            state.layers["hour_angle"] = not state.layers["hour_angle"]
            state.projection_dirty = True
    return True


def run(state: FK6SkyState, width: int, height: int, fixed_time: datetime | None = None) -> None:
    pygame.init()
    pygame.display.set_caption("FK6 Horizon Sky Map")

    #screen = pygame.display.set_mode((width, height))
    screen = pygame.display.set_mode((width, height), pygame.RESIZABLE)

    font, hud_font = pygame.font.SysFont("Microsoft YaHei", 16), pygame.font.SysFont("Microsoft YaHei", 15)
    clock = pygame.time.Clock()
    cached_surface: pygame.Surface | None = None
    next_coordinate_refresh = time.monotonic() + 1.0
    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.VIDEORESIZE:
                width = max(event.w, MIN_WINDOW_WIDTH)
                height = max(event.h, MIN_WINDOW_HEIGHT)
                screen = pygame.display.set_mode((width, height), pygame.RESIZABLE)
                state.projection_dirty = True
                cached_surface = None
                
            running = apply_event(event, state) and running
        now = time.monotonic()
        if fixed_time is None and now >= next_coordinate_refresh:
            state.update_altaz(datetime.now(timezone.utc))
            next_coordinate_refresh = now + 1.0
        if state.projection_dirty or cached_surface is None:
            cached_surface, _ = render_frame(state, width, height, font, hud_font)
        screen.blit(cached_surface, (0, 0))
        draw_star_tooltip(screen, state, font, pygame.mouse.get_pos())
        pygame.display.flip()
        clock.tick(60)
    pygame.quit()


def run_star_map(
    altitude_deg: float | None,
    azimuth_deg: float | None,
    *, longitude_deg: float = 121.55, latitude_deg: float = 29.87,
    magnitude_limit: float | None = None, fixed_time: datetime | None = None,
    export: Path | None = None,
) -> None:
    catalog = ROOT / "data" / "fk6_stars.bin"

    state = FK6SkyState(
        source=catalog,
        latitude_deg=latitude_deg,
        longitude_deg=longitude_deg,
        magnitude_limit=magnitude_limit,
        motion_epoch=fixed_time or datetime.now(timezone.utc),
    )

    if altitude_deg is not None and azimuth_deg is not None:
        state.set_telescope_pointing(altitude_deg, azimuth_deg)

    if export is not None:
        pygame.init()
        font = pygame.font.Font(None, 18)
        surface, _ = render_frame(state, 1200, 900, font, font)
        export = Path(export)
        export.parent.mkdir(parents=True, exist_ok=True)
        pygame.image.save(surface, str(export))
        pygame.quit()
    else:
        run(state, 1200, 900, fixed_time=fixed_time)
