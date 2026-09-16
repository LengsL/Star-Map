"""Download the two FK6 tables from CDS/VizieR and save them as ECSV."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
import ssl
import urllib.request

import certifi
from astropy.table import Table


VIZIER_URL = "https://vizier.cds.unistra.fr/viz-bin/votable?-source={table}&-out.all&-out.max=10000"
TABLES = {
    "I/264/fk6_1": ("fk6_part1.ecsv", 878),
    "I/264/fk6_3": ("fk6_part3.ecsv", 3272),
}


def download_table(table_id: str, destination: Path, expected_rows: int) -> None:
    """Fetch one FK6 VizieR table and write an ECSV file with full metadata."""
    if destination.exists():
        try:
            catalogue = Table.read(destination, format="ascii.ecsv")
            if len(catalogue) == expected_rows:
                print(f"{table_id}: using local {destination} ({len(catalogue)} rows)")
                return
            print(f"{table_id}: local {destination} has {len(catalogue)} rows; redownloading")
        except Exception as error:
            print(f"{table_id}: local {destination} is unreadable ({error}); redownloading")
    request = urllib.request.Request(
        VIZIER_URL.format(table=table_id),
        headers={"User-Agent": "tianguangsuo-fk6-catalogue/1.0"},
    )
    context = ssl.create_default_context(cafile=certifi.where())
    with urllib.request.urlopen(request, context=context, timeout=90) as response:
        payload = response.read()
    catalogue = Table.read(BytesIO(payload), format="votable")
    if len(catalogue) != expected_rows:
        raise RuntimeError(f"{table_id}: expected {expected_rows} rows, received {len(catalogue)}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    catalogue.write(destination, format="ascii.ecsv", overwrite=True)
    print(f"{table_id}: {len(catalogue)} rows -> {destination}")


def main() -> None:
    raw_dir = Path(__file__).resolve().parents[1] / "data" / "raw"
    for table_id, (filename, expected_rows) in TABLES.items():
        download_table(table_id, raw_dir / filename, expected_rows)


if __name__ == "__main__":
    main()
