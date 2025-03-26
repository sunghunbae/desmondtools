import os
import sys
import argparse
import subprocess

from pathlib import Path
from importlib.resources import files

script_path = files('desmondtools.script')

SCHRODINGER = os.environ["SCHRODINGER"]
USER=os.environ["USER"]

schrodinger_hosts = "{}/schrodinger.hosts".format(SCHRODINGER)
schrodinger_run = "{}/run".format(SCHRODINGER)
multisim = "{}/utilities/multisim".format(SCHRODINGER)

# Below python scripts were written with Schrodinger Python API
# and have be wrapped with $SCHRODINGER/run

def batch_cmsinfo() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-cmsinfo.py')] + argv)

def batch_dihedral() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-dihedral.py')] + argv)

def batch_distance() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-distance.py')] + argv)


def batch_ligrmsd() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-ligrmsd.py')] + argv)


def batch_maeinfo() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-maeinfo.py')] + argv)


def batch_mmgbsa() -> None:
    parser = argparse.ArgumentParser(description="trajectory thermal_mmgbsa",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ligand', dest="ligand", default="chain. A", help="ASL for ligand")
    parser.add_argument('--solvent', dest="solvent", default="res. SPC T3P NA CL", help="ASL for solvent")
    parser.add_argument('--start', dest="start", type=int, help="start frame")
    parser.add_argument('--end', dest="end", type=int, help="end frame")
    parser.add_argument('--step', dest="step", type=int, default=1, help="frame step size")
    parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
    args = parser.parse_args()
    for cms_file in args.cms:
        jobname = Path(cms_file).name
        if jobname.startswith("desmond_md_job"):
            jobname = jobname.replace("desmond_md_job", "")
        if jobname.endswith("-out.cms"):
            jobname = jobname.replace("-out.cms", "")
        if not Path(f"{jobname}-prime-out.csv").exists():
            print(f"{cms_file} {args.start} {args.end} --> {jobname}")
            subprocess.run([schrodinger_run,
                            "thermal_mmgbsa.py",
                            "-j", jobname,
                            "-HOST", "localhost:10",
                            "-NJOBS", "2",
                            "-lig_asl", args.ligand,
                            "-atom_asl", args.solvent,
                            "-start_frame", str(args.start),
                            "-end_frame", str(args.end),
                            "-step_size", str(args.step),
                            cms_file,
                            ])
        else:
            print(f"{cms_file} {args.start} {args.end} --> Done/Skip")


def batch_report() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-report.py')] + argv)


def batch_rg() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-rg.py')] + argv)

def batch_rmsx() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-rmsx.py')] + argv)


