"""Interactive FK6 horizon sky map with cached coordinate updates.

The update path is intentionally separated:

* FK6 is read and J2000 proper motion is propagated once at startup;
* local Alt/Az (including the hour-angle grid) is refreshed once per second;
* zoom, constrained pan and magnitude controls merely reproject cached Alt/Az values.
"""

from __future__ import annotations

import argparse
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

# J2000 RA/Dec (degrees) for familiar asterism lines.  Keeping these compact
# coordinates locally makes the FK6 renderer self-contained: it never reads a
# Hipparcos catalogue or a Hipparcos/Stellarium constellation JSON at runtime.
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
    """Catalogue and coordinate caches for the interactive renderer."""

    source: Path
    latitude_deg: float = 29.87
    longitude_deg: float = 121.55
    magnitude_limit: float | None = None
    motion_epoch: datetime | None = None

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
        self.layers = {
            "sun": True, "moon": True, "planets": True, "constellations": True,
            "north_pole": True, "equator": True, "ecliptic": True,
        }
        self.constellation_edges = load_major_constellation_lines()
        self.constellation_altaz: list[tuple[str, np.ndarray, np.ndarray]] = []
        self.observation_time = self.motion_epoch
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
        for hour in (-6, -4, -2, 0, 2, 4, 6):
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
        """Keep the fixed horizon circle fully inside the zoomed sky disc."""
        maximum = self.zoom - 1.0
        distance = math.hypot(*self.pan)
        if distance > maximum and distance > 0.0:
            scale = maximum / distance
            self.pan = (self.pan[0] * scale, self.pan[1] * scale)

    def pan_by(self, dx: float, dy: float) -> None:
        """Move the view, never farther than the current zoom can support."""
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


def draw_horizon_grid(surface: pygame.Surface, font: pygame.font.Font,
                      centre: tuple[int, int], radius: int, labels: LabelPlacer) -> None:
    """Draw the fixed altitude/azimuth grid and physical horizon first."""
    for altitude_deg in (30, 60):
        ring_radius = round(radius * (90 - altitude_deg) / 90)
        pygame.draw.circle(surface, GRID, centre, ring_radius, 1)
        label = font.render(f"Alt {altitude_deg}°", True, GRID_LABEL)
        labels.place(surface, label, (centre[0] + 5, centre[1] - ring_radius - 16))
    for azimuth_deg in range(0, 360, 30):
        x, y = project_altaz(0.0, math.radians(azimuth_deg), centre, radius, 1.0, (0.0, 0.0))
        pygame.draw.line(surface, GRID, centre, (round(float(x)), round(float(y))), 1)
    pygame.draw.circle(surface, HORIZON, centre, radius, 2)
    for azimuth_deg, label in ((0, "N"), (90, "E"), (180, "S"), (270, "W")):
        x, y = project_altaz(0.0, math.radians(azimuth_deg), centre, radius + 19, 1.0, (0.0, 0.0))
        text = font.render(label, True, HORIZON)
        labels.place(surface, text, (round(float(x)), round(float(y))), centered=True)
    pygame.draw.circle(surface, HORIZON, centre, 4, 1)
    zenith = font.render("Zenith", True, HORIZON)
    labels.place(surface, zenith, (centre[0], centre[1] + 18), centered=True)


def draw_hour_angle_lines(surface: pygame.Surface, state: FK6SkyState,
                          font: pygame.font.Font, centre: tuple[int, int], radius: int,
                          labels: LabelPlacer) -> None:
    """Draw cached hour-angle curves after the horizon coordinate grid."""
    for hour, altitude, azimuth in state.hour_angle_lines:
        x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom,
                             (state.pan[0] * radius, state.pan[1] * radius))
        inside = ((altitude >= 0.0) & ((x - centre[0]) ** 2 + (y - centre[1]) ** 2 <= radius ** 2))
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


def _view_pan(state: FK6SkyState, radius: int) -> tuple[float, float]:
    return state.pan[0] * radius, state.pan[1] * radius


def draw_cached_curve(surface: pygame.Surface, altitude: np.ndarray, azimuth: np.ndarray,
                      colour: tuple[int, int, int], state: FK6SkyState,
                      centre: tuple[int, int], radius: int) -> None:
    """Draw the visible parts of a cached Alt/Az curve inside the horizon."""
    x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom, _view_pan(state, radius))
    inside = ((altitude >= 0.0) & ((x - centre[0]) ** 2 + (y - centre[1]) ** 2 <= radius ** 2))
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
    if state.layers["equator"]:
        altitude, azimuth = state.reference_altaz["equator"]
        draw_cached_curve(surface, altitude, azimuth, EQUATOR, state, centre, radius)
    if state.layers["ecliptic"]:
        altitude, azimuth = state.reference_altaz["ecliptic"]
        draw_cached_curve(surface, altitude, azimuth, ECLIPTIC, state, centre, radius)
    if state.layers["constellations"]:
        for _, altitude, azimuth in state.constellation_altaz:
            draw_cached_curve(surface, altitude, azimuth, CONSTELLATION, state, centre, radius)


def draw_stars(surface: pygame.Surface, state: FK6SkyState,
               centre: tuple[int, int], radius: int) -> int:
    """Draw FK6 stars; their appearance is determined only by their own Vmag."""
    x, y = project_altaz(state.altitude_rad, state.azimuth_rad, centre, radius, state.zoom,
                         _view_pan(state, radius))
    magnitudes = state.stars["vmag"]
    inside = ((x - centre[0]) ** 2 + (y - centre[1]) ** 2 <= radius ** 2)
    visible = ((state.altitude_rad >= 0.0) & np.isfinite(magnitudes)
               & (magnitudes <= state.magnitude_limit) & inside)
    for index in np.flatnonzero(visible):
        brightness, dot_radius = star_appearance(float(magnitudes[index]))
        point = (round(float(x[index])), round(float(y[index])))
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


def draw_object_layers(surface: pygame.Surface, state: FK6SkyState,
                       font: pygame.font.Font, centre: tuple[int, int], radius: int,
                       labels: LabelPlacer) -> None:
    """Draw NCP then Sun, Moon and planets on top of the stellar catalogue."""
    if state.layers["north_pole"]:
        altitude, azimuth = state.north_pole_altaz
        if altitude >= 0.0:
            x, y = project_altaz(altitude, azimuth, centre, radius, state.zoom, _view_pan(state, radius))
            point = (round(float(x)), round(float(y)))
            if math.dist(point, centre) <= radius:
                pygame.draw.circle(surface, NORTH_POLE, point, 6, 1)
                pygame.draw.line(surface, NORTH_POLE, (point[0] - 8, point[1]), (point[0] + 8, point[1]), 1)
                pygame.draw.line(surface, NORTH_POLE, (point[0], point[1] - 8), (point[0], point[1] + 8), 1)
                labels.place_near(surface, font.render("NCP", True, NORTH_POLE), point,
                                  ((9, -17), (9, 7), (-38, -17), (-38, 7)))
    for body in state.solar_system:
        if body.name == "Sun" and not state.layers["sun"]:
            continue
        if body.name == "Moon" and not state.layers["moon"]:
            continue
        if body.name not in {"Sun", "Moon"} and not state.layers["planets"]:
            continue
        if body.altitude_rad < 0.0:
            continue
        x, y = project_altaz(body.altitude_rad, body.azimuth_rad, centre, radius, state.zoom,
                             _view_pan(state, radius))
        point = (round(float(x)), round(float(y)))
        if math.dist(point, centre) <= radius:
            pygame.draw.circle(surface, body.colour, point, body.radius)
            labels.place_near(surface, font.render(body.name, True, body.colour), point,
                              ((body.radius + 4, -10), (body.radius + 4, 7),
                               (-body.radius - 50, -10), (-body.radius - 50, 7)))


def star_appearance(magnitude: float) -> tuple[int, int]:
    """Stable Vmag-to-dot mapping, independent of the selected Vmag limit."""
    level = max(0.0, min(1.0, (8.5 - magnitude) / 10.0))
    brightness = round(105 + 150 * (level ** 0.45))
    radius = (4 if magnitude <= 0.0 else 3 if magnitude <= 2.0 else
              2 if magnitude <= 4.0 else 1 if magnitude <= 6.0 else 0)
    return brightness, radius


def render_frame(state: FK6SkyState, width: int, height: int,
                 font: pygame.font.Font, hud_font: pygame.font.Font) -> tuple[pygame.Surface, int]:
    """Reproject only already-cached Alt/Az coordinates onto a new surface."""
    surface = pygame.Surface((width, height))
    surface.fill(BACKGROUND)
    map_left = min(width - 260, SIDEBAR_WIDTH + SIDEBAR_GAP)
    map_width = width - map_left - 18
    centre = (map_left + map_width // 2, height // 2 + 24)
    radius = max(100, min(map_width // 2 - 20, height // 2 - 62))
    labels = LabelPlacer()
    draw_horizon_grid(surface, font, centre, radius, labels)
    draw_hour_angle_lines(surface, state, font, centre, radius, labels)
    draw_reference_layers(surface, state, centre, radius)
    visible = draw_stars(surface, state, centre, radius)
    draw_object_layers(surface, state, font, centre, radius, labels)

    panel = pygame.Rect(16, 16, SIDEBAR_WIDTH - 32, 286)
    panel_surface = pygame.Surface(panel.size, pygame.SRCALPHA)
    pygame.draw.rect(panel_surface, (*HUD_BACKGROUND, HUD_ALPHA), panel_surface.get_rect(), border_radius=8)
    pygame.draw.rect(panel_surface, (*HUD_BORDER, 225), panel_surface.get_rect(), 1, border_radius=8)
    surface.blit(panel_surface, panel)
    details = (
        "FK6 · Horizon / Hour-angle sky map",
        f"Time (UTC): {state.observation_time:%Y-%m-%d %H:%M:%S}",
        f"Location: {state.longitude_deg:.3f}° E, {state.latitude_deg:.3f}° N",
        f"Vmag limit: {state.magnitude_limit:.1f}  (catalog {state.catalog_mag_min:.1f}–{state.catalog_mag_max:.1f})",
        f"Visible: {visible:,} / {len(state.stars):,} FK6 stars",
        f"Zoom ×{state.zoom:.2f}; FOV {180.0 / state.zoom:.1f}°; view {state.pan[0]:+.2f}, {state.pan[1]:+.2f}",
        "Proper motion: startup once; Alt/Az refresh: 1 s",
        f"Coordinate refresh count: {state.coordinate_update_count}",
        "Wheel zoom · arrows / WASD view",
        "[ ] Vmag · 0 reset · Esc quit",
        "1 Sun · 2 Moon · 3 planets",
        "4 Con · 5 NCP · 6 Eq · 7 Ecl",
    )
    for line_number, detail in enumerate(details):
        colour = HUD_TEXT if line_number == 0 else GRID_LABEL
        surface.blit(hud_font.render(detail, True, colour), (panel.x + 12, panel.y + 10 + 21 * line_number))
    state.projection_dirty = False
    return surface, visible


def parse_time(value: str | None) -> datetime | None:
    if value is None:
        return None
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return parsed.replace(tzinfo=timezone.utc) if parsed.tzinfo is None else parsed


def apply_event(event: pygame.event.Event, state: FK6SkyState) -> bool:
    """Apply display-only controls. Return ``False`` when the app should quit."""
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
            state.layers["sun"] = not state.layers["sun"]
            state.projection_dirty = True
        elif event.key == pygame.K_2:
            state.layers["moon"] = not state.layers["moon"]
            state.projection_dirty = True
        elif event.key == pygame.K_3:
            state.layers["planets"] = not state.layers["planets"]
            state.projection_dirty = True
        elif event.key == pygame.K_4:
            state.layers["constellations"] = not state.layers["constellations"]
            state.projection_dirty = True
        elif event.key == pygame.K_5:
            state.layers["north_pole"] = not state.layers["north_pole"]
            state.projection_dirty = True
        elif event.key == pygame.K_6:
            state.layers["equator"] = not state.layers["equator"]
            state.projection_dirty = True
        elif event.key == pygame.K_7:
            state.layers["ecliptic"] = not state.layers["ecliptic"]
            state.projection_dirty = True
    return True


def run(state: FK6SkyState, width: int, height: int, fixed_time: datetime | None = None) -> None:
    pygame.init()
    pygame.display.set_caption("FK6 Horizon Sky Map")
    screen = pygame.display.set_mode((width, height))
    font, hud_font = pygame.font.SysFont("Microsoft YaHei", 16), pygame.font.SysFont("Microsoft YaHei", 15)
    clock = pygame.time.Clock()
    cached_surface: pygame.Surface | None = None
    next_coordinate_refresh = time.monotonic() + 1.0
    running = True
    while running:
        for event in pygame.event.get():
            running = apply_event(event, state) and running
        now = time.monotonic()
        if fixed_time is None and now >= next_coordinate_refresh:
            state.update_altaz(datetime.now(timezone.utc))
            next_coordinate_refresh = now + 1.0
        if state.projection_dirty or cached_surface is None:
            cached_surface, _ = render_frame(state, width, height, font, hud_font)
        screen.blit(cached_surface, (0, 0))
        pygame.display.flip()
        clock.tick(60)
    pygame.quit()


def main() -> None:
    parser = argparse.ArgumentParser(description="Interactive FK6 horizon sky map")
    parser.add_argument("--catalog", type=Path, default=ROOT / "data" / "fk6_stars.bin")
    parser.add_argument("--longitude", type=float, default=121.55, help="east-positive degrees")
    parser.add_argument("--latitude", type=float, default=29.87, help="north-positive degrees")
    parser.add_argument("--magnitude-limit", type=float)
    parser.add_argument("--time", help="ISO 8601 fixed time; primarily for reproducible export")
    parser.add_argument("--export", type=Path, help="write one PNG instead of opening a window")
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()
    fixed_time = parse_time(args.time)
    state = FK6SkyState(args.catalog, args.latitude, args.longitude,
                        args.magnitude_limit, fixed_time or datetime.now(timezone.utc))
    if args.export:
        pygame.init()
        font, hud_font = pygame.font.SysFont("Microsoft YaHei", 16), pygame.font.SysFont("Microsoft YaHei", 15)
        image, visible = render_frame(state, args.width, args.height, font, hud_font)
        args.export.parent.mkdir(parents=True, exist_ok=True)
        pygame.image.save(image, args.export)
        pygame.quit()
        print(f"Saved {args.export} ({visible} visible FK6 stars)")
        return
    run(state, args.width, args.height, fixed_time)


if __name__ == "__main__":
    main()
