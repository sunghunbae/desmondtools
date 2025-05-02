import sys
import argparse

from desmondtools import Maestro


def export():
    parser = argparse.ArgumentParser(description="Convert Maestro file(s) to PDB/mmcif file(s)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--skip-first', dest="skip_first", default=False, action="store_true",
                        help="skip the first entry (i.e. receptor)")
    parser.add_argument('--only-first', dest="only_first", default=False, action="store_true",
                        help="skip the first entry (i.e. receptor)")
    parser.add_argument('--names', dest="names", default=[], nargs="+",
                        help="separately write entries with given name(s)")
    parser.add_argument('--as-complex', dest="as_complex", default=False, action="store_true",
                        help="output as a complex (1 + n-th entries together)")
    parser.add_argument('--separately', dest="separately", default=False, action="store_true",
                        help="output as separate files")
    parser.add_argument('--pdb', dest="pdb", default=False, action="store_true",
                        help="output in PDB format")
    parser.add_argument('--cif', dest="cif", default=False, action="store_true",
                        help="output in mmcif format")
    parser.add_argument('--dict', dest="dict", default=False, action="store_true",
                        help="output in dictionary")
    parser.add_argument('mae', nargs='+', default=[], help='input maestro filename(s)')
    args = parser.parse_args() 
    # args is a Namespace object and can be unpacked into a dictionary using the vars()
    # vars() returns __dict__ attribute of an object.

    if len(args.mae) == 0:
        argparse.print_help()
        sys.exit(0)
    
    for filename in args.mae:
        if args.dict:
            print(Maestro(filename).to_dict())
        else:
            Maestro(filename).export2(**vars(args))