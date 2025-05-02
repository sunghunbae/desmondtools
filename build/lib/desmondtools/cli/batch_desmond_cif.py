import sys
import argparse

from desmondtools import Maestro


def export():
    parser = argparse.ArgumentParser(description="Convert Maestro file(s) to PDB/SDF/mmcif file(s)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--skip-first', dest="skip_first", default=False, action="store_true",
                        help="skip the first entry (i.e. receptor)")
    parser.add_argument('--only-first', dest="only_first", default=False, action="store_true",
                        help="skip the first entry (i.e. receptor)")
    parser.add_argument('--names', dest="names", default=[], nargs="+",
                        help="include given entry name(s)")
    parser.add_argument('--as-complex', dest="as_complex", default=False, action="store_true",
                        help="output as a complex (1 + n-th entries together)")
    parser.add_argument('--separately', dest="separately", default=False, action="store_true",
                        help="output as separate files")
    parser.add_argument('--pdb', dest="pdb", default=False, action="store_true",
                        help="output in PDB format")
    parser.add_argument('--sdf', dest="sdf", default=False, action="store_true",
                        help="output in SDF format")
    parser.add_argument('--mmcif', dest="mmcif", default=False, action="store_true",
                        help="output in mmcif format")
    parser.add_argument('mae', nargs='+', default=[], help='input maestro filename(s)')
    args = parser.parse_args() 
    # args is a Namespace object and can be unpacked into a dictionary using the vars()
    # vars() returns __dict__ attribute of an object.

    if len(args.mae) == 0:
        argparse.print_help()
        sys.exit(0)
    
    for filename in args.mae:
        if args.pdb :
            Maestro(filename).to_pdb(**vars(args))



def batch_cif():
    parser = argparse.ArgumentParser(description="Convert Maestro file(s) to CIF file(s)",
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