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
        Maestro(filename).to_mmcif()


def batch_pdb():
    parser = argparse.ArgumentParser(description="COnvert Maestro file(s) to PDB file(s)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mae', nargs='+', default=[], help='input maestro filename(s)')
    args = parser.parse_args()

    if len(args.mae) == 0:
        argparse.print_help()
        sys.exit(0)
    
    for filename in args.mae:
        Maestro(filename).to_pdb()


def batch_sdf():
    parser = argparse.ArgumentParser(description="COnvert Maestro file(s) to SDF file(s)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mae', nargs='+', default=[], help='input maestro filename(s)')
    args = parser.parse_args()

    if len(args.mae) == 0:
        argparse.print_help()
        sys.exit(0)
    
    for filename in args.mae:
        Maestro(filename).to_sdf()