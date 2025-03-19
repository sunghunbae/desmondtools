import sys
import os
import argparse
import pandas as pd

try:
    from schrodinger import structure
except ImportError:
    print("Schrodinger Python API is required")
    sys.exit(0)

parser = argparse.ArgumentParser(description="Trajectory RMSD and RMSF",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--eaf', dest="eaf", default="", help="event analysis outfile")
parser.add_argument('--prefix', dest="prefix", default="", help="output prefix")
parser.add_argument('cms', help="desmond cms output file")
args = parser.parse_args()

st_reader = structure.StructureReader(args.cms)
st = next(st_reader)

tbl = {}
rmsd_data = {'time_ns':[], 'rmsd':[]}
rmsf_data = {"resSeq":[], "resName":[], "name":[], "atom_index":[], "rmsf":[]}

with open(args.eaf, "r") as f:
    for line in f:
        line = line.strip()
        if "Name" in line:
            name = line.split()[2]
            tbl[name] = None
        
        if "Result" in line:
            data = line.split("=")[1]
            data = data.replace("[","")
            data = data.replace("]","")
            data = data.split()
            tbl[name] = data
        
        if "AtomIndices" in line:
            data = line.split("=")[1]
            data = data.replace("[","")
            data = data.replace("]","")
            data = data.split()
            tbl["AtomIndices"] = data
        
        if "TrajectoryFirstTime" in line:
            TrajectoryFirstTime = float(line.split()[2])
        
        if "TrajectoryInterval_ps" in line:
            TrajectoryInterval_ps = float(line.split()[2])
        
        if "TrajectoryNumFrames" in line:
            TrajectoryNumFrames = int(line.split()[2])
    
    for i, x in enumerate(tbl["RMSD_A"]):
        rmsd_str = f"{float(x):.3f}"
        rmsd_data['time_ns'].append(f"{TrajectoryFirstTime + i*TrajectoryInterval_ps*0.001:.3f}")
        rmsd_data['rmsd'].append(rmsd_str)
    
    for i, x in enumerate(tbl["RMSF_A"]):
        a = st.atom[int(tbl["AtomIndices"][i])]
        rmsf_str = f"{float(x):.3f}"
        rmsf_data['resSeq'].append(a.resnum)
        rmsf_data['resName'].append(a.pdbres.strip())
        rmsf_data['name'].append(a.pdbname.strip())
        rmsf_data['atom_index'].append(a.index)
        rmsf_data['rmsf'].append(rmsf_str)

rmsd_df = pd.DataFrame(rmsd_data)
rmsf_df = pd.DataFrame(rmsf_data).sort_values(by=['resSeq','resName','atom_index'])

rmsd_df.to_csv(f"{args.prefix}-rmsd.csv", index=False)
rmsf_df.to_csv(f"{args.prefix}-rmsf.csv", index=False)