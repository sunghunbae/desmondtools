import sys
import os
import argparse

# entry > molecule > chain > residue > atom

try:
    from schrodinger import structure
    from schrodinger.structutils import analyze
    from schrodinger.application.desmond.packages import traj_util
except ImportError:
    print("Schrodinger Python API is required")
    sys.exit(0)

parser = argparse.ArgumentParser(description="mae information",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--receptor', dest="receptor", default="", help="ASL for receptor")
# parser.add_argument('--ligand', dest="ligand", default="", help="ASL for ligand")
parser.add_argument('mae', nargs="+", help="mae or maegz file(s)")
args = parser.parse_args()

# examples:
#   receptor_asl = '(chain.name "A" AND res.num 833 AND atom.name "SG")'
#   ligand_asl = 'smarts. [N]C(=O)[OH1]'

indent = 4

for mae_file in args.mae:
    for st_index, st in enumerate(structure.StructureReader(mae_file), start=1):
        try:
            title = st.property['s_m_title']
        except:
            title = None
        try:
            entry_id = st.property['s_m_entry_id']
        except:
            entry_id = None
        try:
            entry_name = st.property['s_m_entry_name']
        except:
            entry_name = None

        print(f"(mol. {st_index}) {title}")
        for chain in st.chain:
            print(" " * indent, end='')
            print(f"(chain. {chain.name})")
            residue_numbers = None
            residue_names = {}
            for residue in chain.residue:
                first_atom = list(residue.atom)[0]
                resSeq = first_atom.property['i_m_residue_number']
                resName = first_atom.property['s_m_pdb_residue_name']
                if residue_numbers is None:
                    residue_numbers = [ resSeq ]
                else:
                    residue_numbers.append( resSeq )
                residue_names[resSeq] = resName

            start = None
            last = None
            continous_residues = []
            for resSeq in sorted(residue_numbers):
                if start is None:
                    start = resSeq
                if (last is not None) and (resSeq-last) != 1 :
                    continous_residues.append((start, last))
                    start = resSeq
                last = resSeq
            continous_residues.append((start, last))

            for i, (start, end) in enumerate(continous_residues):
                if i > 0:
                    (prev_start, prev_end) = continous_residues[i-1]
                    print(" " * indent*2, end='')
                    print(f"{prev_end+1:4d} .. {start-1:4d} ({start-prev_end-1}) * missing *")    
                print(" " * indent*2, end='')
                print(f"{residue_names[start]}{start:4d} .. ", end="")
                print(f"{residue_names[end]}{end:4d} ({end-start+1})")
        print()

    # receptor_atoms = []
    # ligands = analyze.find_ligands(st)
    # print(f"number of ligands: {len(ligands)}")
    # num_entry = 0
    # for mol in st:
    #     if not receptor_atoms :
    #         receptor_atoms = analyze.evaluate_asl(mol, args.receptor)
    #         continue
    #     # ligand entries
    #     ligand_atoms = analyze.evaluate_asl(mol, args.ligand)
    #     print(receptor_atoms, ligand_atoms)
    # print("total", num_entry,"molecules")