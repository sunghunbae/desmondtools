import os
import sys
import subprocess

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

def batch_maeinfo() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-maeinfo.py')] + argv)

def batch_report() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-report.py')] + argv)

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

def batch_rg() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-rg.py')] + argv)

def batch_rmsx() -> None:
    argv = sys.argv[1:]
    subprocess.run([schrodinger_run, 
                    script_path.joinpath('batch-desmond-rmsx.py')] + argv)


