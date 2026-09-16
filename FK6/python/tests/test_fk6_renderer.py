from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pygame

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from renderer.projection import (BACKGROUND, CONSTELLATION_NAMES, DECLINATION_GRID_DEGREES,
                                 FK6SkyState, HORIZON, HOUR_ANGLE_GRID_HOURS, LabelPlacer,
                                 MAX_PAN, MAJOR_ASTERISMS, MIN_WINDOW_HEIGHT, MIN_WINDOW_WIDTH,
                                 MoonPhaseInfo, SIDEBAR_GAP, SIDEBAR_WIDTH, TELESCOPE,
                                 TelescopePointing, apply_event,
                                 find_hovered_star, project_altaz, render_frame, star_appearance,
                                 view_geometry)


class FK6RendererCacheTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.start = datetime(2026, 7, 31, 12, tzinfo=timezone.utc)
        cls.state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=cls.start)

    def test_startup_propagates_all_stars_once(self) -> None:
        self.assertEqual(len(self.state.stars), 4150)
        self.assertEqual(self.state.propagation_count, 1)
        self.assertEqual(len(self.state.altitude_rad), 4150)
        self.assertEqual(len(self.state.hour_angle_lines), 24)
        self.assertEqual(tuple(hour for hour, _, _ in self.state.hour_angle_lines), HOUR_ANGLE_GRID_HOURS)

    def test_celestial_layers_are_cached_with_coordinate_refresh(self) -> None:
        self.assertEqual(len(self.state.reference_altaz["equator"][0]), 361)
        self.assertEqual(len(self.state.reference_altaz["ecliptic"][0]), 361)
        for declination_deg in DECLINATION_GRID_DEGREES:
            altitude, azimuth = self.state.reference_altaz[f"dec_{declination_deg:+d}"]
            self.assertEqual(len(altitude), 361)
            self.assertEqual(len(azimuth), 361)
        self.assertGreater(len(self.state.constellation_edges), 0)
        names = {body.name for body in self.state.solar_system}
        self.assertIn("Sun", names)
        self.assertIn("Moon", names)
        self.assertIsInstance(self.state.moon_phase, MoonPhaseInfo)
        self.assertGreaterEqual(self.state.moon_phase.illumination, 0.0)
        self.assertLessEqual(self.state.moon_phase.illumination, 1.0)
        self.assertIsInstance(self.state.moon_phase.waxing, bool)

    def test_every_major_constellation_has_a_display_name(self) -> None:
        self.assertEqual(set(CONSTELLATION_NAMES), set(MAJOR_ASTERISMS))

    def test_controls_only_mark_projection_dirty(self) -> None:
        before = self.state.coordinate_update_count
        self.state.projection_dirty = False
        self.state.zoom_by(1.2)
        self.state.pan_by(0.1, -0.1)
        self.state.adjust_magnitude(-0.25)
        self.assertEqual(self.state.coordinate_update_count, before)
        self.assertEqual(self.state.propagation_count, 1)
        self.assertTrue(self.state.projection_dirty)

    def test_display_layer_keys_are_projection_only(self) -> None:
        state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=self.start)
        before = state.coordinate_update_count
        expected_layers = {
            "constellations": True, "equatorial": True,
            "altaz_grid": True, "hour_angle": True,
        }
        self.assertEqual(state.layers, expected_layers)
        for key, layer_name in (
            (pygame.K_1, "constellations"),
            (pygame.K_2, "equatorial"),
            (pygame.K_3, "altaz_grid"),
            (pygame.K_4, "hour_angle"),
        ):
            state.projection_dirty = False
            self.assertTrue(apply_event(pygame.event.Event(pygame.KEYDOWN, key=key), state))
            self.assertFalse(state.layers[layer_name])
            self.assertTrue(state.projection_dirty)
        self.assertEqual(state.coordinate_update_count, before)
        self.assertEqual(state.propagation_count, 1)

    def test_horizon_frame_remains_when_altaz_grid_is_hidden(self) -> None:
        pygame.font.init()
        state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=self.start)
        state.layers["altaz_grid"] = False
        width, height = 900, 700
        font = pygame.font.Font(None, 18)
        surface, _ = render_frame(state, width, height, font, font)
        map_left = min(width - 260, SIDEBAR_WIDTH + SIDEBAR_GAP)
        map_width = width - map_left - 18
        centre = (map_left + map_width // 2, height // 2 + 24)
        radius = max(100, min(map_width // 2 - 20, height // 2 - 62))
        top_horizon_pixel = surface.get_at((centre[0], centre[1] - radius))[:3]
        self.assertEqual(top_horizon_pixel, HORIZON)

    def test_zenith_marker_follows_altaz_grid_visibility(self) -> None:
        pygame.font.init()
        state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=self.start)
        state.layers.update({
            "altaz_grid": False, "hour_angle": False,
            "equatorial": False, "constellations": False,
        })
        state.magnitude_limit = state.catalog_mag_min
        width, height = 900, 700
        font = pygame.font.Font(None, 18)
        surface, _ = render_frame(state, width, height, font, font)
        map_left = min(width - 260, SIDEBAR_WIDTH + SIDEBAR_GAP)
        map_width = width - map_left - 18
        centre = (map_left + map_width // 2, height // 2 + 24)
        self.assertEqual(surface.get_at(centre)[:3], BACKGROUND)

    def test_pan_moves_the_complete_view_and_has_a_finite_limit(self) -> None:
        self.state.reset_projection()
        self.state.pan_by(0.5, -0.25)
        self.assertEqual(self.state.pan, (0.5, -0.25))
        view_centre, view_radius = view_geometry(self.state, (100, 100), 80)
        zenith_x, zenith_y = project_altaz(np.pi / 2.0, 0.0, (100, 100), 80,
                                           self.state.zoom, (40.0, -20.0))
        self.assertEqual((round(float(zenith_x)), round(float(zenith_y))),
                         (round(view_centre[0]), round(view_centre[1])))
        self.assertEqual(view_radius, 80.0)
        self.state.pan_by(2.0, 2.0)
        self.assertLessEqual(np.hypot(*self.state.pan), MAX_PAN + 1e-12)

    def test_one_second_refresh_recalculates_altaz_without_repropagating(self) -> None:
        before = self.state.coordinate_update_count
        self.state.update_altaz(self.start + timedelta(seconds=1))
        self.assertEqual(self.state.coordinate_update_count, before + 1)
        self.assertEqual(self.state.propagation_count, 1)
        self.assertTrue(np.all(self.state.altitude_rad >= -np.pi / 2.0))
        self.assertTrue(np.all(self.state.altitude_rad <= np.pi / 2.0))

    def test_projection_uses_horizon_coordinates(self) -> None:
        x, y = project_altaz(np.array([np.pi / 2.0, 0.0]), np.array([0.0, 0.0]),
                             (100, 100), 80, 1.0, (0.0, 0.0))
        self.assertEqual((round(x[0]), round(y[0])), (100, 100))
        self.assertEqual((round(x[1]), round(y[1])), (100, 20))

    def test_minimum_resizable_window_can_render(self) -> None:
        pygame.font.init()
        font = pygame.font.Font(None, 18)
        surface, visible = render_frame(self.state, MIN_WINDOW_WIDTH, MIN_WINDOW_HEIGHT, font, font)
        self.assertEqual(surface.get_size(), (MIN_WINDOW_WIDTH, MIN_WINDOW_HEIGHT))
        self.assertGreaterEqual(visible, 0)

    def test_telescope_pointing_loads_from_file_and_draws(self) -> None:
        pygame.font.init()
        with tempfile.TemporaryDirectory() as temporary_directory:
            telescope_file = Path(temporary_directory) / "telescope_pointing.json"
            telescope_file.write_text(json.dumps({"alt_deg": 90.0, "az_deg": 0.0}), encoding="utf-8")
            state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=self.start,
                                telescope_file=telescope_file)
            self.assertEqual(state.telescope_pointing, TelescopePointing(90.0, 0.0))
            font = pygame.font.Font(None, 18)
            surface, _ = render_frame(state, 900, 700, font, font)
            map_left = min(900 - 260, SIDEBAR_WIDTH + SIDEBAR_GAP)
            map_width = 900 - map_left - 18
            centre = (map_left + map_width // 2, 700 // 2 + 24)
            self.assertEqual(surface.get_at(centre)[:3], TELESCOPE)

    def test_hovered_star_uses_cached_screen_positions(self) -> None:
        pygame.font.init()
        font = pygame.font.Font(None, 18)
        _, visible = render_frame(self.state, 900, 700, font, font)
        self.assertGreater(visible, 0)
        self.assertGreater(len(self.state.visible_star_hits), 0)
        first = self.state.visible_star_hits[0]
        hit = find_hovered_star(self.state, (first.x + 1, first.y + 1))
        self.assertIsNotNone(hit)
        self.assertEqual(hit.index, first.index)
        self.assertIsNone(find_hovered_star(self.state, (-1000, -1000)))

    def test_star_appearance_is_independent_of_magnitude_limit(self) -> None:
        before = star_appearance(4.5)
        self.state.adjust_magnitude(-1.0)
        self.assertEqual(before, star_appearance(4.5))
        self.assertGreaterEqual(star_appearance(8.0)[0], 105)
        self.assertLessEqual(star_appearance(-1.0)[0], 255)
        self.assertEqual(star_appearance(0.0)[1], 4)
        self.assertEqual(star_appearance(2.0)[1], 3)
        self.assertEqual(star_appearance(4.0)[1], 2)
        self.assertEqual(star_appearance(6.0)[1], 1)
        self.assertEqual(star_appearance(6.1)[1], 0)

    def test_overlapping_labels_are_not_drawn_twice(self) -> None:
        pygame.font.init()
        surface = pygame.Surface((100, 60))
        label = pygame.font.Font(None, 20).render("Moon", True, (255, 255, 255))
        placer = LabelPlacer()
        self.assertTrue(placer.place(surface, label, (10, 10)))
        self.assertFalse(placer.place(surface, label, (12, 12)))


if __name__ == "__main__":
    unittest.main()
