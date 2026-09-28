from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
MISSING = "—"
@dataclass(frozen=True)
class StarInfo:
    id: str
    fk6: str = ""
    hip: str = ""
    hd: str = ""
    hr: str = ""
    name: str = ""
    spectral_type: str = ""
    bv: str = ""
    parallax_mas: float | None = None
    distance_pc: float | None = None
    constellation: str = ""
    @property
    def display_name(self):
        return self.name.strip() or (f"HIP {self.hip}" if self.hip else f"FK6 {self.fk6}")
    def value(self, field):
        value=getattr(self, field)
        return MISSING if value in (None, "") else str(value)
def load_metadata(source: Path, expected_count: int):
    if not source.exists(): return [StarInfo(f"FK6 {i+1}", fk6=str(i+1)) for i in range(expected_count)]
    payload=json.loads(source.read_text(encoding="utf-8")); rows=payload.get("stars", payload)
    if len(rows)!=expected_count: raise ValueError(f"{source}: metadata count mismatch")
    return [StarInfo(**row) for row in rows]
