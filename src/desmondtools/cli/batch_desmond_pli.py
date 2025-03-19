import argparse
import sys
import pandas as pd

from desmondtools import Interaction


def batch_pli():
    teststring ="""Keywords = [
    {RMSD = {
        ASL = "((mol. 1 and backbone) and not (atom.ele H) and not (mol. 2))"
        Frame = 0
        Panel = pl_interact_survey
        Result = [0.0 1.161 1.286 1.331 1.176 1.195 ]
        SelectionType = Backbone
        Tab = pl_rmsd_tab
        Type = ASL
        Unit = Angstrom
        }
    }
    {RMSD = {
        ASL = "(mol. 1 and sidechain and not (mol. 2))"
        Frame = 0
        Panel = pl_interact_survey
        Result = [0.0 1.161 1.286 1.331 1.176 1.195 ]
        SelectionType = "Side chains"
        Tab = pl_rmsd_tab
        Type = ASL
        Unit = Angstrom
        }
    }
    ]"""
    # result = expr.parse_string(teststring)
    # d = result.as_dict()
    # traverse_dict(d)
    # dot = DotMap(d)
    # try:
    #     assert dot.Keywords[1].RMSD.SelectionType == '"Side chains"'
    #     print("ok")
    # except:
    #     print("error")

    parser = argparse.ArgumentParser(description="Average Protein-Ligand Interactions",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--out', dest='out', default='mean-PLIF', help="output basename")
    parser.add_argument('eaf', nargs='+', default=[], help='input -out.eaf filename(s)')
    args = parser.parse_args()

    if len(args.eaf) == 0:
        argparse.print_help()
        sys.exit(0)

    total_num_frames = 0
    data = {}
    for filename in args.eaf:
        event_analysis = Interaction(filename, verbose=False)
        total_num_frames += event_analysis.num_frames
        
        print(f"{filename}  {event_analysis.num_frames} frames")
        
        for resSeq in event_analysis.HBond:
            if resSeq in data:
                data[resSeq]['hbond'] += event_analysis.HBond[resSeq]['count']
            else:
                data[resSeq] = {'resName': event_analysis.HBond[resSeq]['resName'], 
                                'hbond': event_analysis.HBond[resSeq]['count'],
                                'hydrophobic': 0,
                                'polar': 0,
                                'waterbridge': 0,
                                }
        for resSeq in event_analysis.Hydrophobic:
            if resSeq in data:
                data[resSeq]['hydrophobic'] += event_analysis.Hydrophobic[resSeq]['count']
            else:
                data[resSeq] = {'resName': event_analysis.Hydrophobic[resSeq]['resName'], 
                                'hbond': 0,
                                'hydrophobic': event_analysis.Hydrophobic[resSeq]['count'],
                                'polar': 0,
                                'waterbridge': 0,
                                }
        for resSeq in event_analysis.Polar:
            if resSeq in data:
                data[resSeq]['polar'] += event_analysis.Polar[resSeq]['count']
            else:
                data[resSeq] = {'resName': event_analysis.Polar[resSeq]['resName'], 
                                'hbond': 0,
                                'hydrophobic': 0,
                                'polar': event_analysis.Polar[resSeq]['count'],
                                'waterbridge': 0,
                                }
        for resSeq in event_analysis.WaterBridge:
            if resSeq in data:
                data[resSeq]['waterbridge'] += event_analysis.WaterBridge[resSeq]['count']
            else:
                data[resSeq] = {'resName': event_analysis.WaterBridge[resSeq]['resName'],
                                'hbond': 0,
                                'hydrophobic' : 0,
                                'polar': 0,
                                'waterbridge': event_analysis.WaterBridge[resSeq]['count'],
                                } 

    csvdata = {'resid':[], 
            'resSeq':[], 
            'resName':[], 
            'hbond':[], 
            'hydrophobic':[],
            'polar': [],
            'waterbridge': [],
            }
    
    for resSeq in sorted(data):
        csvdata['resSeq'].append(resSeq)
        csvdata['resName'].append(data[resSeq]['resName'])
        csvdata['resid'].append(f"{data[resSeq]['resName']}_{resSeq}")
        csvdata['hbond'].append(float(data[resSeq]['hbond'])/total_num_frames)
        csvdata['hydrophobic'].append(float(data[resSeq]['hydrophobic'])/total_num_frames)
        csvdata['polar'].append(float(data[resSeq]['polar'])/total_num_frames)
        csvdata['waterbridge'].append(float(data[resSeq]['waterbridge'])/total_num_frames)

    df = pd.DataFrame(csvdata)
    df.to_csv(args.out + '.csv', index=False, float_format='%.4f')
    g = df.loc[:, ~df.columns.isin(['resSeq','resName'])].plot.bar(
        x="resid", 
        stacked=True,
        title="Protein-Ligand Interactions", 
        xlabel='Residue', 
        ylabel='Fraction of MD trajectory', 
        figsize=(8,3), 
        fontsize=8,
        )
    fig = g.get_figure()
    fig.savefig(args.out + '.pdf', bbox_inches="tight", pad_inches=0.2, dpi=150)