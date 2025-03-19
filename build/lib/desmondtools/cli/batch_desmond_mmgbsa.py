import argparse
import subprocess
import os

from pathlib import Path


def batch_mmgbsa() -> None:
    SCHRODINGER = os.environ["SCHRODINGER"]
    schrodinger_run = "{}/run".format(SCHRODINGER)

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