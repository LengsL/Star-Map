from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from catalog.fk6_reader import iter_stars, read_header


class FK6CatalogueTest(unittest.TestCase):
    def test_merged_binary_has_all_records(self) -> None:
        source = ROOT / "data" / "fk6_stars.bin"
        self.assertEqual(read_header(source), 4150)
        self.assertEqual(sum(1 for _ in iter_stars(source)), 4150)


if __name__ == "__main__":
    unittest.main()
