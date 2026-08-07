"""Load the 88 Western IAU constellation stick figures keyed by HIP number."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class Constellation:
    abbreviation: str
    name: str
    edges: tuple[tuple[int, int], ...]


def load_constellations(source: Path) -> tuple[Constellation, ...]:
    """Load Stellarium Modern IAU JSON and turn polylines into HIP edges."""
    document = json.loads(source.read_text(encoding="utf-8"))
    constellations: list[Constellation] = []
    for item in document["constellations"]:
        edges: list[tuple[int, int]] = []
        for polyline in item.get("lines", []):
            edges.extend(zip(polyline, polyline[1:]))
        name = item.get("common_name", {}).get("native") or item["id"]
        constellations.append(Constellation(item["id"].rsplit(" ", 1)[-1], name, tuple(edges)))
    if len(constellations) != 88:
        raise ValueError(f"expected 88 IAU constellations, received {len(constellations)}")
    return tuple(constellations)
