
import argparse
from pathlib import Path
from renderer.projection import parse_time
from renderer.projection import run_star_map

def main():

    parser = argparse.ArgumentParser(description="Interactive FK6 horizon sky map")

    parser.add_argument("--altitude", type=float, default=None,)
    parser.add_argument("--azimuth", type=float, default=None)
    parser.add_argument("--longitude", type=float, default=121.55)
    parser.add_argument("--latitude", type=float, default=29.87)
    parser.add_argument("--magnitude-limit", type=float, default=None)
    parser.add_argument("--time", type=parse_time, default=None)
    parser.add_argument("--export", type=Path, default=None)

    args = parser.parse_args()

    if (args.altitude is None) != (args.azimuth is None):
        parser.error("--altitude and --azimuth must be supplied together")
    run_star_map(args.altitude, args.azimuth, longitude_deg=args.longitude,
                 latitude_deg=args.latitude, magnitude_limit=args.magnitude_limit,
                 fixed_time=args.time,
                 export=args.export)




if __name__ == "__main__":
    main()
