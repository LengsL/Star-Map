
import argparse
from pathlib import Path
from renderer.projection import (ObservatoryLocation, load_observatory_config,
                                 parse_time, run_star_map)

def main():

    parser = argparse.ArgumentParser(description="Interactive FK6 horizon sky map")

    parser.add_argument("--altitude", type=float, default=None,)
    parser.add_argument("--azimuth", type=float, default=None)
    parser.add_argument("--longitude", type=float, default=None,
                        help="east-positive longitude in degrees; overrides observatory config")
    parser.add_argument("--latitude", type=float, default=None,
                        help="north-positive latitude in degrees; overrides observatory config")
    parser.add_argument("--height", type=float, default=None,
                        help="observatory height above WGS84 ellipsoid in metres")
    parser.add_argument("--observatory", type=Path, default=None,
                        help="path to observatory JSON (longitude_deg, latitude_deg, height_m)")
    parser.add_argument("--goto-output", type=Path, default=None,
                        help="write selected-object GOTO commands to this JSON file")
    parser.add_argument("--magnitude-limit", type=float, default=None)
    parser.add_argument("--time", type=parse_time, default=None)
    parser.add_argument("--export", type=Path, default=None)

    args = parser.parse_args()

    if (args.altitude is None) != (args.azimuth is None):
        parser.error("--altitude and --azimuth must be supplied together")
    configured = load_observatory_config(args.observatory)
    observatory = ObservatoryLocation(
        longitude_deg=configured.longitude_deg if args.longitude is None else args.longitude,
        latitude_deg=configured.latitude_deg if args.latitude is None else args.latitude,
        height_m=configured.height_m if args.height is None else args.height,
    )
    run_star_map(args.altitude, args.azimuth, observatory=observatory,
                 magnitude_limit=args.magnitude_limit,
                 fixed_time=args.time,
                 export=args.export,
                 goto_output=args.goto_output)




if __name__ == "__main__":
    main()
