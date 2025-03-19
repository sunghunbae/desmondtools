import argparse
import subprocess
import os

from desmondtools import Multisim

from importlib.resources import files

import desmondtools
script_path = files('desmondtools.script')


SCHRODINGER = os.environ["SCHRODINGER"]

USER=os.environ["USER"]

schrodinger_hosts = "{}/schrodinger.hosts".format(SCHRODINGER)
schrodinger_run = "{}/run".format(SCHRODINGER)

multisim = "{}/utilities/multisim".format(SCHRODINGER)

def batch_rmsx() -> None:
    parser = argparse.ArgumentParser(description="trajectory RMSD and RMSF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--protein', dest="protein", default="(mol. 1) and not atom.ele H", help="ASL for protein")
    parser.add_argument('--ligand', dest="ligand", default="(mol. 2) and not atom.ele H", help="ASL for ligand")
    parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
    args = parser.parse_args()
    event_analysis_infile  = "temporary-event-in.eaf"
    event_analysis_outfile = "temporary-event-out.eaf"
    events = Multisim(template="desmond-rmsx-analysis.eaf")
    events.dot.Keywords[0].RMSD.ASL = args.protein
    events.dot.Keywords[1].RMSF.ASL = args.ligand
    events.dot.Keywords[1].RMSF.FitBy = args.protein
    events.dot.ProteinASL = args.protein
    events.dot.LigandASL = args.ligand
    for cms_file in args.cms:
        print(f"analyzing {cms_file} ...")
        prefix = os.path.basename(cms_file).replace('-out.cms','')
        trj = cms_file.replace("-out.cms", "_trj")
        events.dot.Trajectory = cms_file
        with open(event_analysis_infile, "w") as f:
            events.write(f)
        # requires Schrodinger Python API
        subprocess.run([schrodinger_run, 
                        "analyze_simulation.py",
                        cms_file, 
                        trj, 
                        event_analysis_outfile, 
                        event_analysis_infile],
                        stdout=subprocess.DEVNULL, 
                        stderr=subprocess.STDOUT)
        if os.path.exists(event_analysis_infile):
            os.remove(event_analysis_infile)
        # requires Schrodinger Python API
        subprocess.run([schrodinger_run, 
                        script_path.joinpath('batch-desmond-rmsx.py'),
                        cms_file, 
                        event_analysis_outfile, 
                        prefix])
        if os.path.exists(event_analysis_outfile):
            os.remove(event_analysis_outfile)