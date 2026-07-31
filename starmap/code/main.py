"""Interactive black-sky map backed by the supplied JPL DE files."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import sys

#import pygame

from astronomy import degree, gnomonic, julian_date, radec
from catalog import CONSTELLATIONS, STARS, Star
from jpl_ephemeris import DateOutOfRange, EphemerisError, JPLEphemeris

BODY_NAMES = {0: ("Mercury", -1.0), 1: ("Venus", -4.0), 3: ("Mars", -2.0),
              4: ("Jupiter", -2.7), 5: ("Saturn", 0.5), 6: ("Uranus", 5.7),
              7: ("Neptune", 7.8), 8: ("Pluto", 14.5), 9: ("Moon", -12.0),
              10: ("Sun", -26.7)}


def source_candidates(source_dir: Path) -> list[Path]:
    return sorted([*source_dir.glob("*.200"), *source_dir.glob("*.406")])


def choose_ephemeris(paths: list[Path], jd: float) -> JPLEphemeris | None:
    for path in paths:
        ephemeris = JPLEphemeris(path)
        if ephemeris.header.start_jd <= jd <= ephemeris.header.end_jd:
            return ephemeris
        ephemeris.close()
    return None


class SkyMap:
    def __init__(self, source_dir: Path):
        pygame.init()
        self.screen = pygame.display.set_mode((1280, 800), pygame.RESIZABLE)
        pygame.display.set_caption("天光所 · JPL 星图")
        self.font = pygame.font.SysFont("Microsoft YaHei", 16)
        self.small = pygame.font.SysFont("Microsoft YaHei", 13)
        self.clock = pygame.time.Clock()
        self.files = source_candidates(source_dir)
        self.now = datetime.now(timezone.utc)
        self.jd = julian_date(self.now)
        self.ephemeris = choose_ephemeris(self.files, self.jd)
        self.center_ra, self.center_dec = math.radians(90), math.radians(10)
        self.scale = 420.0
        self.magnitude_limit = 3.5
        self.drag_anchor: tuple[int, int] | None = None
        self.hovered: tuple[str, str, float, float, float] | None = None
        self.status = self._status()

    def _status(self) -> str:
        if self.ephemeris is None:
            listed = ", ".join(p.name for p in self.files) or "未找到 .200/.406 文件"
            return f"当前日期不在可用历表范围内（{listed}）"
        h = self.ephemeris.header
        return f"DE{h.de_number} · JD {self.jd:.2f} · {self.now:%Y-%m-%d %H:%M UTC}"

    def reset_time(self) -> None:
        if self.ephemeris:
            self.ephemeris.close()
        self.now, self.jd = datetime.now(timezone.utc), julian_date(datetime.now(timezone.utc))
        self.ephemeris = choose_ephemeris(self.files, self.jd)
        self.status = self._status()

    def screen_point(self, ra: float, dec: float) -> tuple[int, int] | None:
        point = gnomonic(ra, dec, self.center_ra, self.center_dec)
        if point is None:
            return None
        w, h = self.screen.get_size()
        return int(w / 2 + point[0] * self.scale), int(h / 2 + point[1] * self.scale)

    def body_objects(self) -> list[tuple[str, str, float, float, float]]:
        if not self.ephemeris:
            return []
        objects = []
        for body, (name, magnitude) in BODY_NAMES.items():
            try:
                # One light-time iteration is sufficient for this visual display.
                vector = self.ephemeris.position(self.jd, body, 2)
                distance = math.sqrt(sum(v * v for v in vector))
                vector = self.ephemeris.position(self.jd - distance * 0.0057755, body, 2)
                ra, dec = radec(vector)
                objects.append((name, "JPL DE ephemeris", magnitude, ra, dec))
            except (DateOutOfRange, EphemerisError):
                continue
        return objects

    def draw(self) -> None:
        self.screen.fill((0, 0, 0))
        self.hovered = None
        by_name = {star.name: star for star in STARS}
        # Constellations are deliberately behind stars and planets.
        for _, names in CONSTELLATIONS:
            points = [self.screen_point(by_name[name].ra, by_name[name].dec) for name in names]
            visible = [point for point in points if point]
            if len(visible) > 1:
                pygame.draw.lines(self.screen, (185, 32, 32), False, visible, 1)
        mouse = pygame.mouse.get_pos()
        for star in STARS:
            if star.magnitude > self.magnitude_limit:
                continue
            point = self.screen_point(star.ra, star.dec)
            if not point:
                continue
            radius = max(1, int(4.8 - max(star.magnitude, -1.5)))
            pygame.draw.circle(self.screen, (225, 232, 255), point, radius)
            if math.dist(point, mouse) < max(9, radius + 5):
                self.hovered = (star.name, f"Star · {star.constellation}", star.magnitude, star.ra, star.dec)
        for item in self.body_objects():
            name, kind, magnitude, ra, dec = item
            point = self.screen_point(ra, dec)
            if not point:
                continue
            colour = (255, 217, 110) if name == "Sun" else (130, 210, 255)
            pygame.draw.circle(self.screen, colour, point, 6)
            self.screen.blit(self.small.render(name, True, colour), (point[0] + 8, point[1] - 8))
            if math.dist(point, mouse) < 12:
                self.hovered = item
        self._label(self.status, 12, 10)
        self._label(f"显示星等 ≤ {self.magnitude_limit:.1f}   滚轮缩放 · 拖动平移 · [ / ] 调星等 · T 当前时刻", 12, 32)
        if self.hovered:
            name, kind, mag, ra, dec = self.hovered
            text = f"{name} | {kind} | 星等 {mag:.2f} | RA {degree(ra)/15:.2f}h · Dec {degree(dec):.2f}°"
            surface = self.font.render(text, True, (255, 255, 255))
            rect = surface.get_rect(topleft=(mouse[0] + 14, mouse[1] + 14))
            pygame.draw.rect(self.screen, (28, 28, 38), rect.inflate(10, 8), border_radius=4)
            self.screen.blit(surface, rect)
        pygame.display.flip()

    def _label(self, text: str, x: int, y: int) -> None:
        self.screen.blit(self.small.render(text, True, (185, 185, 185)), (x, y))

    def run(self) -> None:
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
                elif event.type == pygame.MOUSEWHEEL:
                    self.scale = min(4000, max(60, self.scale * (1.18 ** event.y)))
                elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                    self.drag_anchor = event.pos
                elif event.type == pygame.MOUSEBUTTONUP and event.button == 1:
                    self.drag_anchor = None
                elif event.type == pygame.MOUSEMOTION and self.drag_anchor:
                    dx, dy = event.pos[0] - self.drag_anchor[0], event.pos[1] - self.drag_anchor[1]
                    self.center_ra = (self.center_ra - dx / self.scale) % (2 * math.pi)
                    self.center_dec = min(math.radians(85), max(math.radians(-85), self.center_dec + dy / self.scale))
                    self.drag_anchor = event.pos
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_LEFTBRACKET:
                        self.magnitude_limit = max(-1.5, self.magnitude_limit - 0.5)
                    elif event.key == pygame.K_RIGHTBRACKET:
                        self.magnitude_limit = min(8.0, self.magnitude_limit + 0.5)
                    elif event.key == pygame.K_r:
                        self.center_ra, self.center_dec, self.scale = math.radians(90), math.radians(10), 420.0
                    elif event.key == pygame.K_t:
                        self.reset_time()
            self.draw()
            self.clock.tick(30)
        if self.ephemeris:
            self.ephemeris.close()
        pygame.quit()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="天光所 JPL 星图")
    parser.add_argument("--source-dir", type=Path,
                        default=Path(os.getenv("JPL_SOURCE_DIR", r"E:\UNNC\天光所\1\sourece")),
                        help="包含 lnxm*.200 和 lnxm*.406 的目录")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not args.source_dir.is_dir():
        sys.exit(f"找不到 JPL 数据目录：{args.source_dir}")
    SkyMap(args.source_dir).run()
