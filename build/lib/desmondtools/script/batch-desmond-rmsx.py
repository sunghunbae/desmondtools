import sys
import os
import pandas as pd

try:
    from schrodinger import structure
except ImportError:
    print("schrodinger python API is required")
    sys.exit(0)

(cms_file, event_analysis_outfile, prefix) = (sys.argv[1], sys.argv[2], sys.argv[3])

st_reader = structure.StructureReader(cms_file)
st = next(st_reader)

tbl = {}
rmsd_data = {'time_ns':[], 'rmsd':[]}
rmsf_data = {"resSeq":[], "resName":[], "name":[], "atom_index":[], "rmsf":[]}

with open(event_analysis_outfile, "r") as f:
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
rmsd_df.to_csv(f"{prefix}-rmsd.csv", index=False)
rmsf_df.to_csv(f"{prefix}-rmsf.csv", index=False)