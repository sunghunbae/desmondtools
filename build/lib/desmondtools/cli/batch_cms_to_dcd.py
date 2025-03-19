import argparse
import subprocess
import desmondtools

from pathlib import Path
from importlib.resources import files


script_path = files('desmondtools.script')


def cms_to_dcd() -> None:
    parser = argparse.ArgumentParser(description="convert desmond to dcd",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
    args = parser.parse_args()
    for cms_file in args.cms:
        outdcd = Path(cms_file).stem + "-nowat.dcd"
        if not Path(outdcd).exists():
            subprocess.run(['vmd', 
                            '-display', 'text', 
                            '-e', script_path.joinpath('vmd_convert_desmond_to_dcd.tcl'),
                            '-args', cms_file,
                            ])