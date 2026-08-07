from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys
import unittest

import numpy as np
import pygame

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from renderer.projection import FK6SkyState, LabelPlacer, project_altaz, star_appearance


class FK6RendererCacheTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.start = datetime(2026, 7, 31, 12, tzinfo=timezone.utc)
        cls.state = FK6SkyState(ROOT / "data" / "fk6_stars.bin", motion_epoch=cls.start)

    def test_startup_propagates_all_stars_once(self) -> None:
        self.assertEqual(len(self.state.stars), 4150)
        self.assertEqual(self.state.propagation_count, 1)
        self.assertEqual(len(self.state.altitude_rad), 4150)
        self.assertEqual(len(self.state.hour_angle_lines), 7)

    def test_celestial_layers_are_cached_with_coordinate_refresh(self) -> None:
        self.assertEqual(len(self.state.reference_altaz["equator"][0]), 361)
        self.assertEqual(len(self.state.reference_altaz["ecliptic"][0]), 361)
        self.assertGreater(len(self.state.constellation_edges), 0)
        names = {body.name for body in self.state.solar_system}
        self.assertIn("Sun", names)
        self.assertIn("Moon", names)

    def test_controls_only_mark_projection_dirty(self) -> None:
        before = self.state.coordinate_update_count
        self.state.projection_dirty = False
        self.state.zoom_by(1.2)
        self.state.pan_by(0.1, -0.1)
        self.state.adjust_magnitude(-0.25)
        self.assertEqual(self.state.coordinate_update_count, before)
        self.assertEqual(self.state.propagation_count, 1)
        self.assertTrue(self.state.projection_dirty)

    def test_pan_is_limited_by_the_current_zoom(self) -> None:
        self.state.reset_projection()
        self.state.pan_by(1.0, 1.0)
        self.assertEqual(self.state.pan, (0.0, 0.0))
        self.state.zoom_by(1.5)
        self.state.pan_by(2.0, 2.0)
        self.assertLessEqual(np.hypot(*self.state.pan), self.state.zoom - 1.0 + 1e-12)

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
