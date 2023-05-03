from .koffi import ImageMetadata, skybot_search_frame, jpl_search_frame
from astropy.table import QTable
from astropy import units as u
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("filename", type=str, help="The fits filts to open")
    parser.add_argument("--query", type=str, default="JPL", help="The service to query, JPL or SkyBot")
    parser.add_argument("--format", type=str, default="QTable", help="The output format to use")

    args = parser.parse_args()

    image = ImageMetadata(args.filename)

    if args.query == "JPL":
        search = jpl_search_frame
    elif args.query == "SkyBot":
        search = skybot_search_frame

    objects = search(image)
    names = []
    ras = []
    decs = []
    xs = []
    ys = []
    for name, coord in objects:
        names.append(name)
        ras.append(coord.ra.degree * u.degree)
        decs.append(coord.dec.degree * u.degree)
        x, y = image.wcs.world_to_pixel(coord)
        xs.append(x)
        ys.append(y)

    table = QTable(
        data=[names, ras, decs, xs, ys], 
        names=["Object Name", "RA", "Dec", "x", "y"],
    )
    if args.format == "QTable":
        print(table)
    else:
        table.write(sys.stdout, format=args.format)


if __name__ == "__main__":
    main()
