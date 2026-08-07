"""Observer-facing circular sky map from Hipparcos ``stars.bin``.

For a time, longitude and latitude, convert RA/Dec to Alt/Az, exclude objects
below the horizon, then project the visible hemisphere into a circular screen.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
from pathlib import Path
import sys

import pygame

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from astronomy.astronomy import (julian_date, radec_to_altaz, topocentric_vector,
                                 vector_to_radec)  # noqa: E402
from catalog.constellations import Constellation, load_constellations  # noqa: E402
from catalog.star import Star, read_stars_bin  # noqa: E402
from jpl_ephemeris import EphemerisError, JPLEphemeris  # noqa: E402


BACKGROUND = (0, 0, 0)
GRID = (30, 36, 48)
GRID_LABEL = (128, 140, 160)
HORIZON = (100, 145, 190)
HUD_BACKGROUND = (8, 12, 20)
HUD_TEXT = (225, 232, 245)
CONSTELLATION = (190, 50, 65)
EQUATOR = (70, 135, 230)
ECLIPTIC = (230, 175, 65)
SUN = (255, 214, 95)
MOON = (205, 214, 230)


def project_altaz(altitude_rad: float, azimuth_rad: float,
                  centre: tuple[int, int], radius: int, zoom: float = 1.0,
                  pan: tuple[float, float] = (0.0, 0.0)) -> tuple[int, int]:
    """Zenithal equidistant map: zenith at centre and horizon at edge."""
    distance = radius * zoom * (math.pi / 2.0 - altitude_rad) / (math.pi / 2.0)
    return (round(centre[0] + pan[0] + distance * math.sin(azimuth_rad)),
            round(centre[1] + pan[1] - distance * math.cos(azimuth_rad)))


def appearance(magnitude: float, magnitude_limit: float) -> tuple[int, int]:
    """Map Hipparcos magnitude to grayscale brightness and dot radius."""
    fraction = max(0.0, min(1.0, (magnitude_limit - magnitude) / (magnitude_limit + 1.5)))
    return round(65 + 190 * math.sqrt(fraction)), max(1, min(6, round(1 + 5 * fraction)))


def point_in_view(point: tuple[int, int], centre: tuple[int, int], radius: int) -> bool:
    return math.dist(point, centre) <= radius


def sky_point(ra_rad: float, dec_rad: float, jd_utc: float, latitude_rad: float,
              longitude_rad: float, centre: tuple[int, int], map_radius: int,
              zoom: float, pan: tuple[float, float]) -> tuple[int, int] | None:
    """Project an RA/Dec coordinate only if it is above the local horizon."""
    altitude, azimuth = radec_to_altaz(ra_rad, dec_rad, jd_utc, latitude_rad, longitude_rad)
    if altitude < 0.0:
        return None
    point = project_altaz(altitude, azimuth, centre, map_radius, zoom, pan)
    return point if point_in_view(point, centre, map_radius) else None


def draw_curve(surface: pygame.Surface, coordinates: list[tuple[float, float]],
               colour: tuple[int, int, int], jd_utc: float, latitude_rad: float,
               longitude_rad: float, centre: tuple[int, int], map_radius: int,
               zoom: float, pan: tuple[float, float]) -> None:
    """Draw visible continuous portions of a sampled celestial great-circle."""
    segment: list[tuple[int, int]] = []
    for ra_rad, dec_rad in coordinates:
        point = sky_point(ra_rad, dec_rad, jd_utc, latitude_rad, longitude_rad,
                          centre, map_radius, zoom, pan)
        if point is None:
            if len(segment) > 1:
                pygame.draw.aalines(surface, colour, False, segment)
            segment = []
        else:
            segment.append(point)
    if len(segment) > 1:
        pygame.draw.aalines(surface, colour, False, segment)


def ecliptic_coordinates(samples: int = 361) -> list[tuple[float, float]]:
    """J2000 ecliptic, expressed as equatorial RA/Dec coordinates."""
    obliquity = math.radians(23.4392911)
    coordinates = []
    for index in range(samples):
        longitude = 2.0 * math.pi * index / (samples - 1)
        ra_rad = math.atan2(math.sin(longitude) * math.cos(obliquity), math.cos(longitude)) % (2.0 * math.pi)
        dec_rad = math.asin(math.sin(longitude) * math.sin(obliquity))
        coordinates.append((ra_rad, dec_rad))
    return coordinates


def draw_constellations(surface: pygame.Surface, constellations: tuple[Constellation, ...],
                        stars_by_hip: dict[int, Star], jd_utc: float, latitude_rad: float,
                        longitude_rad: float, centre: tuple[int, int], map_radius: int,
                        zoom: float, pan: tuple[float, float]) -> None:
    for constellation in constellations:
        for first_hip, second_hip in constellation.edges:
            first, second = stars_by_hip.get(first_hip), stars_by_hip.get(second_hip)
            if first is None or second is None:
                continue
            first_point = sky_point(first.ra_rad, first.dec_rad, jd_utc, latitude_rad,
                                    longitude_rad, centre, map_radius, zoom, pan)
            second_point = sky_point(second.ra_rad, second.dec_rad, jd_utc, latitude_rad,
                                     longitude_rad, centre, map_radius, zoom, pan)
            if first_point and second_point:
                pygame.draw.aaline(surface, CONSTELLATION, first_point, second_point)


def open_ephemeris(source_dir: Path, jd_utc: float) -> JPLEphemeris | None:
    """Open the first supplied JPL DE file whose coverage includes the date."""
    for path in sorted([*source_dir.glob("*.200"), *source_dir.glob("*.406")]):
        try:
            candidate = JPLEphemeris(path)
            if candidate.header.start_jd <= jd_utc <= candidate.header.end_jd:
                return candidate
            candidate.close()
        except EphemerisError:
            continue
    return None


def draw_solar_body(surface: pygame.Surface, ephemeris: JPLEphemeris, body: int,
                    label: str, colour: tuple[int, int, int], radius: int, jd_utc: float,
                    latitude_rad: float, longitude_rad: float, centre: tuple[int, int],
                    map_radius: int, zoom: float, pan: tuple[float, float],
                    font: pygame.font.Font) -> None:
    try:
        geocentric = ephemeris.position(jd_utc, body, 2)
        ra_rad, dec_rad = vector_to_radec(topocentric_vector(geocentric, jd_utc,
                                                               latitude_rad, longitude_rad))
        point = sky_point(ra_rad, dec_rad, jd_utc, latitude_rad, longitude_rad,
                          centre, map_radius, zoom, pan)
    except EphemerisError:
        return
    if point:
        pygame.draw.circle(surface, colour, point, radius)
        text = font.render(label, True, colour)
        surface.blit(text, (point[0] + radius + 4, point[1] - radius - 4))


def render(stars: list[Star], width: int, height: int, magnitude_limit: float,
           jd_utc: float, latitude_rad: float, longitude_rad: float,
           observation_time: datetime, latitude_deg: float, longitude_deg: float,
           font: pygame.font.Font, hud_font: pygame.font.Font, zoom: float = 1.0,
           pan: tuple[float, float] = (0.0, 0.0),
           layers: dict[str, bool] | None = None,
           constellations: tuple[Constellation, ...] = (),
           stars_by_hip: dict[int, Star] | None = None,
           ephemeris: JPLEphemeris | None = None) -> tuple[pygame.Surface, int]:
    """Render the upper hemisphere; stars at Alt < 0 are deliberately omitted."""
    surface = pygame.Surface((width, height))
    surface.fill(BACKGROUND)
    centre = (width // 2, height // 2 + 28)
    map_radius = min(width, height) // 2 - 82
    layers = layers or {}
    stars_by_hip = stars_by_hip or {}
    # Altitude rings and azimuth spokes are the screen coordinate grid.
    for altitude_deg in (30, 60):
        ring_radius = round(map_radius * (90 - altitude_deg) / 90)
        pygame.draw.circle(surface, GRID, centre, ring_radius, 1)
        surface.blit(font.render(f"Alt {altitude_deg}°", True, GRID_LABEL),
                     (centre[0] + 5, centre[1] - ring_radius - 14))
    for azimuth_deg in range(0, 360, 30):
        point = project_altaz(0.0, math.radians(azimuth_deg), centre, map_radius)
        pygame.draw.line(surface, GRID, centre, point, 1)
        if azimuth_deg not in (0, 90, 180, 270):
            text = font.render(f"{azimuth_deg}°", True, GRID_LABEL)
            surface.blit(text, text.get_rect(center=point))

    # Celestial reference layers are drawn behind the star catalogue.
    if layers.get("equator", True):
        draw_curve(surface, [(2.0 * math.pi * i / 360.0, 0.0) for i in range(361)],
                   EQUATOR, jd_utc, latitude_rad, longitude_rad, centre, map_radius, zoom, pan)
    if layers.get("ecliptic", True):
        draw_curve(surface, ecliptic_coordinates(), ECLIPTIC, jd_utc, latitude_rad,
                   longitude_rad, centre, map_radius, zoom, pan)
    if layers.get("constellations", True):
        draw_constellations(surface, constellations, stars_by_hip, jd_utc, latitude_rad,
                            longitude_rad, centre, map_radius, zoom, pan)

    visible = 0
    for star in stars:
        if star.magnitude > magnitude_limit:
            continue
        altitude, azimuth = radec_to_altaz(star.ra_rad, star.dec_rad, jd_utc,
                                            latitude_rad, longitude_rad)
        if altitude < 0.0:
            continue
        brightness, dot_radius = appearance(star.magnitude, magnitude_limit)
        point = project_altaz(altitude, azimuth, centre, map_radius, zoom, pan)
        # The physical horizon circle is also the circular viewport boundary.
        if math.dist(point, centre) <= map_radius:
            pygame.draw.circle(surface, (brightness, brightness, brightness), point, dot_radius)
            visible += 1

    if layers.get("sun", True) and ephemeris:
        draw_solar_body(surface, ephemeris, 10, "太阳  Sun", SUN, 7, jd_utc,
                        latitude_rad, longitude_rad, centre, map_radius, zoom, pan, hud_font)
    if layers.get("moon", True) and ephemeris:
        draw_solar_body(surface, ephemeris, 9, "月球  Moon", MOON, 6, jd_utc,
                        latitude_rad, longitude_rad, centre, map_radius, zoom, pan, hud_font)
    if layers.get("north_pole", True):
        pole = sky_point(0.0, math.pi / 2.0, jd_utc, latitude_rad, longitude_rad,
                         centre, map_radius, zoom, pan)
        if pole:
            pygame.draw.circle(surface, EQUATOR, pole, 6, 1)
            pygame.draw.line(surface, EQUATOR, (pole[0] - 8, pole[1]), (pole[0] + 8, pole[1]), 1)
            pygame.draw.line(surface, EQUATOR, (pole[0], pole[1] - 8), (pole[0], pole[1] + 8), 1)
            surface.blit(font.render("北天极", True, EQUATOR), (pole[0] + 9, pole[1] - 15))

    # Fixed foreground labels keep the orientation and horizon readable.
    pygame.draw.circle(surface, HORIZON, centre, map_radius, 2)
    for azimuth_deg, label in ((0, "北  N"), (90, "东  E"),
                               (180, "南  S"), (270, "西  W")):
        point = project_altaz(0.0, math.radians(azimuth_deg), centre, map_radius + 16)
        text = hud_font.render(label, True, HORIZON)
        surface.blit(text, text.get_rect(center=point))
    pygame.draw.circle(surface, HORIZON, centre, 4, 1)
    pygame.draw.line(surface, HORIZON, (centre[0] - 8, centre[1]), (centre[0] + 8, centre[1]), 1)
    pygame.draw.line(surface, HORIZON, (centre[0], centre[1] - 8), (centre[0], centre[1] + 8), 1)
    zenith = hud_font.render("天顶  Zenith", True, HORIZON)
    surface.blit(zenith, zenith.get_rect(center=(centre[0], centre[1] + 18)))

    # Persistent basic sky-map interface.
    panel = pygame.Rect(14, 14, 430, 205)
    pygame.draw.rect(surface, HUD_BACKGROUND, panel, border_radius=7)
    pygame.draw.rect(surface, HORIZON, panel, 1, border_radius=7)
    utc_time = observation_time.astimezone(timezone.utc)
    details = (
        "天光所 · 地平星图",
        f"当前时间  {utc_time:%Y-%m-%d %H:%M:%S UTC}",
        f"当前经纬度  东经 {longitude_deg:.3f}°   北纬 {latitude_deg:.3f}°",
        f"星等上限  Hpmag ≤ {magnitude_limit:.1f}",
        f"视野缩放  ×{zoom:.2f}   视场约 {180.0 / zoom:.1f}°",
        f"当前视野内  {visible:,} / {len(stars):,} 颗",
        "滚轮缩放 · 方向键/WASD 移动 · 0 复位 · [ ] 星等",
    )
    details += (
        "[1] Sun  [2] Moon  [3] Constellations",
        "[4] NCP  [5] Equator  [6] Ecliptic",
    )
    for index, detail in enumerate(details):
        colour = HUD_TEXT if index == 0 else GRID_LABEL
        surface.blit(hud_font.render(detail, True, colour), (27, 24 + index * 19))
    return surface, visible


def parse_time(value: str | None) -> datetime:
    """Parse ISO 8601 and treat a timezone-less value as UTC."""
    if value is None:
        return datetime.now(timezone.utc)
    moment = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return moment.replace(tzinfo=timezone.utc) if moment.tzinfo is None else moment.astimezone(timezone.utc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Hipparcos RA/Dec to Alt/Az circular sky map")
    parser.add_argument("--source", type=Path, default=CODE_ROOT / "data" / "stars.bin")
    parser.add_argument("--constellation-source", type=Path,
                        default=CODE_ROOT / "data" / "modern_iau.json")
    parser.add_argument("--jpl-source-dir", type=Path,
                        default=CODE_ROOT.parent / "sourece",
                        help="directory containing lnxm*.200 / lnxm*.406")
    parser.add_argument("--time", type=str,
                        help="ISO 8601 observation time; default: current UTC")
    parser.add_argument("--longitude", type=float, default=121.55,
                        help="observer longitude in degrees; east positive")
    parser.add_argument("--latitude", type=float, default=29.87,
                        help="observer latitude in degrees; north positive")
    parser.add_argument("--magnitude-limit", type=float, default=7.0)
    parser.add_argument("--zoom", type=float, default=1.0,
                        help="initial view zoom, from 1 to 8")
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=800)
    parser.add_argument("--export", type=Path, help="save a PNG and exit")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.source.is_file():
        raise SystemExit(f"Star catalogue not found: {args.source}; run catalog/star.py first")
    if args.width < 320 or args.height < 320:
        raise SystemExit("Canvas must be at least 320 by 320")
    if not -90.0 <= args.latitude <= 90.0:
        raise SystemExit("Latitude must be between -90 and 90 degrees")

    observation_time = parse_time(args.time)
    jd_utc = julian_date(observation_time)
    latitude_rad = math.radians(args.latitude)
    longitude_rad = math.radians(args.longitude)
    stars = list(read_stars_bin(args.source))
    stars_by_hip = {star.hip: star for star in stars}
    constellations = load_constellations(args.constellation_source) if args.constellation_source.is_file() else ()
    ephemeris = open_ephemeris(args.jpl_source_dir, jd_utc) if args.jpl_source_dir.is_dir() else None
    layers = {
        "sun": True,
        "moon": True,
        "constellations": True,
        "north_pole": True,
        "equator": True,
        "ecliptic": True,
    }

    pygame.init()
    grid_font = pygame.font.SysFont("Microsoft YaHei", 14)
    hud_font = pygame.font.SysFont("Microsoft YaHei", 16)
    magnitude_limit = args.magnitude_limit
    zoom = min(8.0, max(1.0, args.zoom))
    pan = (0.0, 0.0)
    sky, visible = render(stars, args.width, args.height, magnitude_limit,
                          jd_utc, latitude_rad, longitude_rad, observation_time,
                          args.latitude, args.longitude, grid_font, hud_font, zoom, pan,
                          layers, constellations, stars_by_hip, ephemeris)
    if args.export:
        args.export.parent.mkdir(parents=True, exist_ok=True)
        pygame.image.save(sky, args.export)
        pygame.quit()
        if ephemeris:
            ephemeris.close()
        print(f"Saved {args.export}: {visible}/{len(stars)} stars in the circular view")
        return

    screen = pygame.display.set_mode((args.width, args.height))
    pygame.display.set_caption("天光所 · 地平星图")
    dirty, running = True, True
    live_time = args.time is None
    last_second = observation_time.replace(microsecond=0)
    clock = pygame.time.Clock()
    while running:
        elapsed = clock.tick(60) / 1000.0
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.MOUSEWHEEL:
                zoom = min(8.0, max(1.0, zoom * (1.15 ** event.y)))
                dirty = True
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_LEFTBRACKET:
                    magnitude_limit = max(-1.5, magnitude_limit - 0.5)
                    dirty = True
                elif event.key == pygame.K_RIGHTBRACKET:
                    magnitude_limit = min(13.0, magnitude_limit + 0.5)
                    dirty = True
                elif event.key == pygame.K_t:
                    live_time = True
                    dirty = True
                elif event.key == pygame.K_0:
                    zoom, pan = 1.0, (0.0, 0.0)
                    dirty = True
                elif event.key == pygame.K_1:
                    layers["sun"] = not layers["sun"]
                    dirty = True
                elif event.key == pygame.K_2:
                    layers["moon"] = not layers["moon"]
                    dirty = True
                elif event.key == pygame.K_3:
                    layers["constellations"] = not layers["constellations"]
                    dirty = True
                elif event.key == pygame.K_4:
                    layers["north_pole"] = not layers["north_pole"]
                    dirty = True
                elif event.key == pygame.K_5:
                    layers["equator"] = not layers["equator"]
                    dirty = True
                elif event.key == pygame.K_6:
                    layers["ecliptic"] = not layers["ecliptic"]
                    dirty = True
        # Keyboard input pans the star field inside the circular viewport.
        pressed = pygame.key.get_pressed()
        move_speed = 260.0 * elapsed
        pan_x, pan_y = pan
        if pressed[pygame.K_LEFT] or pressed[pygame.K_a]:
            pan_x += move_speed
        if pressed[pygame.K_RIGHT] or pressed[pygame.K_d]:
            pan_x -= move_speed
        if pressed[pygame.K_UP] or pressed[pygame.K_w]:
            pan_y += move_speed
        if pressed[pygame.K_DOWN] or pressed[pygame.K_s]:
            pan_y -= move_speed
        if (pan_x, pan_y) != pan:
            # Keep navigation finite; zooming out or pressing 0 restores the full sky.
            maximum_pan = (min(args.width, args.height) // 2 - 82) * (zoom - 0.75)
            pan = (max(-maximum_pan, min(maximum_pan, pan_x)),
                   max(-maximum_pan, min(maximum_pan, pan_y)))
            dirty = True
        if live_time:
            latest_time = datetime.now(timezone.utc)
            if latest_time.replace(microsecond=0) != last_second:
                observation_time = latest_time
                jd_utc = julian_date(observation_time)
                last_second = observation_time.replace(microsecond=0)
                dirty = True
        if dirty:
            sky, visible = render(stars, args.width, args.height, magnitude_limit,
                                  jd_utc, latitude_rad, longitude_rad, observation_time,
                                  args.latitude, args.longitude, grid_font, hud_font, zoom, pan,
                                  layers, constellations, stars_by_hip, ephemeris)
            screen.blit(sky, (0, 0))
            pygame.display.flip()
            dirty = False
    pygame.quit()
    if ephemeris:
        ephemeris.close()


if __name__ == "__main__":
    main()
