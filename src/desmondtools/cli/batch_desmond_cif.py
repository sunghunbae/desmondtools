import sys
import argparse

from desmondtools import Maestro


def batch_cif():
    parser = argparse.ArgumentParser(description="COnvert Maestro file(s) to CIF file(s)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mae', nargs='+', default=[], help='input maestro filename(s)')
    args = parser.parse_args()

    if len(args.mae) == 0:
        argparse.print_help()
        sys.exit(0)
    
    for filename in args.mae:
        u = Maestro(filename)
        u.convert_to_mmcif()