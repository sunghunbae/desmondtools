import sys
import os
import argparse

try:
    from schrodinger import structure
    from schrodinger.structutils import analyze
    from schrodinger.application.desmond.packages import traj_util
except ImportError:
    print("schrodinger python API is required")
    sys.exit(0)

parser = argparse.ArgumentParser(description="trajectory information",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--asl', dest="asl", default="", help="ASL")
parser.add_argument('cms', nargs="+", help="desmond cms output file(s)")
args = parser.parse_args()
   
for cms_file in args.cms:
    # id
    name = os.path.basename(cms_file).replace("-out.cms","").replace("desmond_md_job_","")

    # read .cms and its associated trajectory file
    (msys_model, cms_model, traj) = traj_util.read_cms_and_traj(cms_file)
    
    # number of frames
    num_frames = len(traj)

    # formal charge
    st = structure.StructureReader.read(cms_file)
    num_atoms = 0
    total_charge = 0
    for chain in st.chain:
        for res in chain.residue:
                for a in res.atom:
                    total_charge += int(a.property["i_m_formal_charge"])
                    num_atoms += 1

    print(f"{cms_file} frames {num_frames:4d} atoms {num_atoms} charge {total_charge}")

    if args.asl:
        atom_indice = analyze.evaluate_asl(st, args.asl)
        for i in atom_indice:
            print(i,
                st.atom[i].property["s_m_chain_name"],
                st.atom[i].property['i_m_residue_number'],
                st.atom[i].property["s_m_pdb_residue_name"],
                st.atom[i].property['s_m_pdb_atom_name'],
                st.atom[i].property["r_m_x_coord"],
                st.atom[i].property["r_m_y_coord"],
                st.atom[i].property["r_m_z_coord"],
                )