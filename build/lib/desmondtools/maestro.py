import gzip
import sys
import re
import random
import string
import shlex
import copy
import io

from pdbecif.mmcif_io import CifFileWriter

"""
References:
    G. van Ginkel, et al., 
    PDBeCIF: an open-source mmCIF/CIF parsing and processing package. 
    BMC Bioinformatics 22, 383 (2021).

    - https://github.com/PDBeurope/pdbecif
    - https://pdbeurope.github.io/pdbecif/

"""

class Maestro:
    """ export maestro file in mmcif format """

    def __init__(self, filename:str, max_two_entries:bool=True) -> None:
        """ initialize """

        if filename.endswith('.maegz'):
            self.prefix = filename.replace(".maegz","")
            self.maesto_file = gzip.open(filename, "rt")
        elif filename.endswith('.mae.gz'):
            self.prefix = filename.replace(".mae.gz","")
            self.maesto_file = gzip.open(filename, "rt")
        elif filename.endswith('.mae'):
            self.prefix = filename.replace(".mae","")
            self.maesto_file = open(filename, "rt")
        else:
            print(".mae, .mae.gz, or .maegz are expected")
            sys.exit(0)
        
        # standard DNA/RNA/Protein
        self.std_residues = [
            "ADE", "GUA", "CYT", "URA", "THY",
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ]

        # conversion to standard residue name
        self.std_rename = {
            "URI" : "URA",
            "CYX" : "CYS", # disulfide-bonded
            "CYM" : "CYS", # deprotonated (- charge) and/or bound to metal atoms
            "HID" : "HIS", # protonated at delta position
            "HIE" : "HIS", # protonated at epsilon
            "HIP" : "HIS", # protonated at both delta and epsilon
        } 

        self.atomic_symbol = {
             1 : "H",   2 : "He",  3 : "Li",  4 : "Be",  5 : "B",   6 : "C",   7 : "N",   8 : "O",   9 : "F",  10 : "Ne",
            11 : "Na", 12 : "Mg", 13 : "Al", 14 : "Si", 15 : "P",  16 : "S",  17 : "Cl", 18 : "Ar", 19 : "K",  20 : "Ca",
            21 : "Sc", 22 : "Ti", 23 : "V",  24 : "Cr", 25 : "Mn", 26 : "Fe", 27 : "Co", 28 : "Ni", 29 : "Cu", 30 : "Zn",
            31 : "Ga", 32 : "Ge", 33 : "As", 34 : "Se", 35 : "Br", 36 : "Kr", 37 : "Rb", 38 : "Sr", 39 : "Y",  40 : "Zr",
            41 : "Nb", 42 : "Mo", 43 : "Tc", 44 : "Ru", 45 : "Rh", 46 : "Pd", 47 : "Ag", 48 : "Cd", 49 : "In", 50 : "Sn",
            51 : "Sb", 52 : "Te", 53 : "I",  54 : "Xe", 55 : "Cs", 56 : "Ba", 57 : "La", 58 : "Ce", 59 : "Pr", 60 : "Nd",
            61 : "Pm", 62 : "Sm", 63 : "Eu", 64 : "Gd", 65 : "Tb", 66 : "Dy", 67 : "Ho", 68 : "Er", 69 : "Tm", 70 : "Yb",
            71 : "Lu", 72 : "Hf", 73 : "Ta", 74 : "W",  75 : "Re", 76 : "Os", 77 : "Ir", 78 : "Pt", 79 : "Au", 80 : "Hg",
            81 : "Tl", 82 : "Pb", 83 : "Bi", 84 : "Po", 85 : "At", 86 : "Rn", 87 : "Fr", 88 : "Ra", 89 : "Ac", 90 : "Th",
            91 : "Pa", 92 : "U",  93 : "Np", 94 : "Pu", 95 : "Am", 96 : "Cm", 97 : "Bk", 98 : "Cf", 99 : "Es", 100: "Fm",
            101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds",
            111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
        }

        self.max_two_entries = max_two_entries
        self.count = re.compile(r'\[(\d+)\]')
        self.title = ""
        self.serial = 0
        self.usedId = {}
        self.mmcif = {
            "_entry": None, 
            "_chem_comp_bond": None, 
            "_atom_site": None,
            } # mmcif-like object
        self.iamap = {}
        self.output = io.StringIO()


    def append_entry(self, title):
        """ start a new entry """

        if self.mmcif["_entry"]:
            if type(self.mmcif["_entry"]["id"]) == list:
                # you got here as the third and later entry
                if self.max_two_entries:
                    # write out last two entries
                    self.write_mmcif()
                    # return to the state in which only the first entry was defined
                    self.mmcif = copy.deepcopy(self.mmcif_first)
                    self.usedId = copy.deepcopy(self.mmcif_first_usedId)
                    self.serial = self.mmcif_first_serial
                    # regardless of its appearance, current entry now is regarded as the second entry
                    self.mmcif["_entry"]["id"] = [ self.mmcif["_entry"]["id"][0], title ]
                else:
                    # otherwise continue to add entry
                    self.mmcif["_entry"]["id"].append(title)
            else:
                # you got here as the second entry
                # save the first entry before adding the second entry
                self.mmcif_first = copy.deepcopy(self.mmcif)
                self.mmcif_first_usedId = copy.deepcopy(self.usedId)
                self.mmcif_first_serial = self.serial
                # change to a list type
                self.mmcif["_entry"]["id"] = [ self.mmcif["_entry"]["id"], title ]
        else:
            # regular string type
            self.mmcif["_entry"] = { "id" : title }


    def append_atom_site(self, 
                        group_PDB, id, 
                        type_symbol, 
                        label_atom_id, 
                        label_alt_id, 
                        label_comp_id,
                        label_asym_id, 
                        label_entity_id, 
                        label_seq_id, 
                        pdbx_PDB_ins_code, 
                        Cartn_x, 
                        Cartn_y, 
                        Cartn_z,
                        occupancy, 
                        B_iso_or_equiv, 
                        pdbx_formal_charge, 
                        auth_seq_id, 
                        auth_comp_id, 
                        auth_asym_id, 
                        auth_atom_id, 
                        pdbx_PDB_model_num):
        
        """ add mmcif _atom_site """

        if self.mmcif["_atom_site"] :
            self.mmcif["_atom_site"]["group_PDB"].append(group_PDB)
            self.mmcif["_atom_site"]["id"].append(id)
            self.mmcif["_atom_site"]["type_symbol"].append(type_symbol)
            self.mmcif["_atom_site"]["label_atom_id"].append(label_atom_id)
            self.mmcif["_atom_site"]["label_alt_id"].append(label_alt_id)
            self.mmcif["_atom_site"]["label_comp_id"].append(label_comp_id)
            self.mmcif["_atom_site"]["label_asym_id"].append(label_asym_id)
            self.mmcif["_atom_site"]["label_entity_id"].append(label_entity_id)
            self.mmcif["_atom_site"]["label_seq_id"].append(label_seq_id)
            self.mmcif["_atom_site"]["pdbx_PDB_ins_code"].append(pdbx_PDB_ins_code)
            self.mmcif["_atom_site"]["Cartn_x"].append(Cartn_x)
            self.mmcif["_atom_site"]["Cartn_y"].append(Cartn_y)
            self.mmcif["_atom_site"]["Cartn_z"].append(Cartn_z)
            self.mmcif["_atom_site"]["occupancy"].append(occupancy)
            self.mmcif["_atom_site"]["B_iso_or_equiv"].append(B_iso_or_equiv)
            self.mmcif["_atom_site"]["pdbx_formal_charge"].append(pdbx_formal_charge)
            self.mmcif["_atom_site"]["auth_seq_id"].append(auth_seq_id)
            self.mmcif["_atom_site"]["auth_comp_id"].append(auth_comp_id)
            self.mmcif["_atom_site"]["auth_asym_id"].append(auth_asym_id)
            self.mmcif["_atom_site"]["auth_atom_id"].append(auth_atom_id)
            self.mmcif["_atom_site"]["pdbx_PDB_model_num"].append(pdbx_PDB_model_num)
        else:
            self.mmcif["_atom_site"] = {
                "group_PDB" : [ group_PDB ],
                "id" : [ id ],
                "type_symbol" : [ type_symbol ],
                "label_atom_id" : [ label_atom_id ],
                "label_alt_id" : [ label_alt_id ],
                "label_comp_id" : [ label_comp_id ],
                "label_asym_id" : [ label_asym_id ],
                "label_entity_id" : [ label_entity_id ],
                "label_seq_id" : [ label_seq_id ],
                "pdbx_PDB_ins_code" : [ pdbx_PDB_ins_code ],
                "Cartn_x" : [ Cartn_x ],
                "Cartn_y" : [ Cartn_y ],
                "Cartn_z" : [ Cartn_z ],
                "occupancy" : [ occupancy ],
                "B_iso_or_equiv" : [ B_iso_or_equiv ],
                "pdbx_formal_charge" : [ pdbx_formal_charge ],
                "auth_seq_id" : [ auth_seq_id ],
                "auth_comp_id" : [ auth_comp_id ],
                "auth_asym_id" : [ auth_asym_id ],
                "auth_atom_id" : [ auth_atom_id ],
                "pdbx_PDB_model_num" : [ pdbx_PDB_model_num ],
            }
    

    def append_chem_comp_bond(self, chem_comp_id, i, j, k):
        """ add chem_comp_bond """
        p = self.mmcif["_atom_site"]["label_atom_id"][self.mmcif["_atom_site"]["id"].index(i)]
        q = self.mmcif["_atom_site"]["label_atom_id"][self.mmcif["_atom_site"]["id"].index(j)]
        if self.mmcif["_chem_comp_bond"]:
            self.mmcif["_chem_comp_bond"]["comp_id"].append(chem_comp_id)
            self.mmcif["_chem_comp_bond"]["atom_id_1"].append(p)
            self.mmcif["_chem_comp_bond"]["atom_id_2"].append(q)
            self.mmcif["_chem_comp_bond"]["value_order"].append(k)
        else:
            self.mmcif["_chem_comp_bond"] = {
                "comp_id" : [ chem_comp_id ],
                "atom_id_1" : [ p ],
                "atom_id_2" : [ q ],
                "value_order" : [ k ],
            }


    def new_chem_comp(self):
        """ start a new chem_comp entry """
        if not type(self.mmcif["_entry"]["id"]) == list or len(self.mmcif["_entry"]["id"]) <= 2:
            self.chem_comp_id = "LIG"
        else:
            self.chem_comp_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
        self.chem_comp_chainId = [ c for c in string.ascii_uppercase if not (c in self.usedId)][0]
        self.chem_comp_resSeq = 100
        self.usedId[self.chem_comp_chainId] = []
        return self.chem_comp_chainId, self.chem_comp_resSeq, self.chem_comp_id
        

    def new_chem_comp_atom(self, obj, element):
        """ create and return unique atom name """
        if "s_m_pdb_atom_name" in obj:
            atom_name = obj["s_m_pdb_atom_name"].strip()
        else:
            atom_name = "{}{}".format(element,1)
        if not (atom_name in self.usedId[self.chem_comp_chainId]):
            self.usedId[self.chem_comp_chainId].append(atom_name)
            return atom_name
        else:
            for i in range(1,10000):
                new_name = "{}{}".format(element,i)
                if not (new_name in self.usedId[self.chem_comp_chainId]):
                    self.usedId[self.chem_comp_chainId].append(new_name)
                    return new_name
    

    def value_or_default(self, obj, k, default):
        """ return value or default """
        v = default
        if k in obj:
            if type(obj[k]) == str:
                v = obj[k].strip()
                if not v:
                    v = default
        return v


    def convert_to_mmcif(self):
        token = None
        entity_id = 0
        for line in self.maesto_file:
            line = line.strip()
            
            if not line or line.startswith("#"): 
                continue
            
            # f_m_ct block
            if line.startswith("f_m_ct {"):
                token = 'f_m_ct_key'
                natoms = 0
                nbonds = 0
                entity_id += 1
                k = []
                continue
            if token == 'f_m_ct_key':
                if line.startswith(":::"): # end of f_m_ct key block
                    n = len(k)
                    token = 'f_m_ct_val'
                    v = []
                    continue
                k.append(line)
            if token == 'f_m_ct_val':
                v.append(line)
                if len(v) == n:
                    data = dict(zip(k,v))
                    if "s_lp_Variant" in data:
                        self.title = data['s_lp_Variant'].replace('"','').replace("'","")
                    else:
                        self.title = data['s_m_title'].replace('"','').replace("'","")
                    self.append_entry(self.title)
                    token = None
                    
            
            # m_atom block
            if line.startswith("m_atom") :
                natoms = int(self.count.findall(line)[0])
                token = 'm_atom_key'
                k = ['i_atom_index'] # 1st column is atom_index
                continue
            if token == 'm_atom_key':
                if line.startswith(":::"):
                    token = 'm_atom_val'
                    nv = 0
                    m_atom = {}
                    m_atom_resName = {}
                    continue
                k.append(line)
            if token == 'm_atom_val':
                v = ['' if x=="<>" else x for x in shlex.split(line)]
                data = dict(zip(k,v))

                # build hierarchical structure of chain/ resSeq/ resName
                chainId = self.value_or_default(data,"s_m_chain_name", "A")
                resSeq = int(data["i_m_residue_number"])
                resName = self.value_or_default(data, "s_m_pdb_residue_name", "?")
                if resName in self.std_rename:
                    resName = self.std_rename[resName]
                if not chainId in m_atom:
                    m_atom[chainId] = {}
                if not resSeq in m_atom[chainId]:
                    m_atom[chainId][resSeq] = {}
                if not resName in m_atom[chainId][resSeq]:
                    m_atom[chainId][resSeq][resName] = []
                m_atom[chainId][resSeq][resName].append(data)
                m_atom_resName[int(data["i_atom_index"])] = resName

                nv += 1
                if nv == natoms: # end of m_atom block
                    token = None
                    
            """ m_bond block """
            if line.startswith("m_bond") :
                nbonds = int(self.count.findall(line)[0])
                token = 'm_bond_key'
                k = ['i_bond_index']
                continue
            if token == 'm_bond_key':
                if line.startswith(":::"): # end of m_bond key block
                    token = 'm_bond_val'
                    nv = 0
                    m_bond = {}
                    continue
                k.append(line)
            if token == 'm_bond_val':
                v = ['' if x=="<>" else x for x in shlex.split(line)]
                data = dict(zip(k,v))

                # build bond connectivity
                p = int(data["i_m_from"]) 
                q = int(data["i_m_to"]) 
                r = int(data["i_m_order"])
                if p in m_bond:
                    m_bond[p][q] = r
                else:
                    m_bond[p] = {q:r}
                if q in m_bond:
                    m_bond[q][p] = r
                else:
                    m_bond[q] = {p:r}
        
                nv += 1
                if nv == nbonds: # m_bond block ends
                    token = None

                    # create mmcif _atom_site
                    for chainId_ in sorted(m_atom):
                        for resSeq_ in sorted(m_atom[chainId_]):
                            for resName_ in sorted(m_atom[chainId_][resSeq_]):
                                if resName_ in self.std_residues:
                                    group_PDB = "ATOM"
                                    chainId, resSeq, resName = chainId_, resSeq_, resName_
                                else:
                                    group_PDB = "HETATM"
                                    # check if new_chem_comp is necessary
                                    this_residue = [ int(d["i_atom_index"]) for d in m_atom[chainId_][resSeq_][resName_] ]
                                    need_new_chem_comp = True
                                    for p_ in this_residue:
                                        for q_, bond_order in m_bond[p_].items():
                                            if  (not (q_ in this_residue)) and \
                                                (not (m_atom_resName[q_] in self.std_residues)) and \
                                                (q_ in self.iamap):
                                                # do not create another chem_comp
                                                need_new_chem_comp = False
                                                q, chainId, resSeq, resName = self.iamap[q_]
                                                break
                                    if need_new_chem_comp:
                                        chainId, resSeq, resName = self.new_chem_comp()
                                
                                # go through all atoms within the same residue
                                for d in m_atom[chainId_][resSeq_][resName_]:
                                    self.serial += 1
                                    self.iamap[int(d["i_atom_index"])] = (self.serial, chainId, resSeq, resName)
                                    element = self.atomic_symbol[int(d["i_m_atomic_number"])]
                                    if group_PDB == "ATOM":
                                        name = self.value_or_default(d,"s_m_pdb_atom_name", "?")
                                    else:
                                        name = self.new_chem_comp_atom(d, element) 
                                    x = float(d["r_m_x_coord"])
                                    y = float(d["r_m_y_coord"])
                                    z = float(d["r_m_z_coord"])
                                    occupancy = float(self.value_or_default(d,"r_m_pdb_occupancy",1))
                                    bfactor = float(self.value_or_default(d,"r_m_pdb_tfactor", 0))
                                    formal_charge = int(self.value_or_default(d,"i_m_formal_charge", 0))
                                    self.append_atom_site(
                                        group_PDB = group_PDB, 
                                        id = self.serial,
                                        type_symbol = element,
                                        label_atom_id = name,
                                        label_alt_id = ".",
                                        label_comp_id = resName,
                                        label_asym_id = chainId,
                                        label_entity_id = entity_id,
                                        label_seq_id = resSeq,
                                        pdbx_PDB_ins_code = "?",
                                        Cartn_x = x,
                                        Cartn_y = y,
                                        Cartn_z = z,
                                        occupancy = occupancy,
                                        B_iso_or_equiv = bfactor,
                                        pdbx_formal_charge = formal_charge,
                                        auth_seq_id = resSeq,
                                        auth_comp_id = resName,
                                        auth_asym_id = chainId,
                                        auth_atom_id = name,
                                        pdbx_PDB_model_num = 1,
                                    )

                    # create mmcif _chem_comp_bond
                    for chainId_ in sorted(m_atom):
                        for resSeq_ in sorted(m_atom[chainId_]):
                            for resName_ in sorted(m_atom[chainId_][resSeq_]):
                                if not (resName_ in self.std_residues): # HETATM
                                    for d in m_atom[chainId_][resSeq_][resName_]:
                                        p_ = int(d["i_atom_index"])
                                        p, p_chainId, p_resSeq, p_resName = self.iamap[p_]
                                        if p_ in m_bond:
                                            for q_, bond_order in m_bond[p_].items():
                                                q, q_chainId, q_resSeq, q_resName = self.iamap[q_]
                                                self.append_chem_comp_bond(p_resName, p, q, bond_order)
                    
        self.maesto_file.close()
        self.write_mmcif()
        print(self.output.getvalue())
        self.output.close()


    def write_mmcif(self, in_memory:bool=True):
        """ write to .cif file """

        if type(self.mmcif["_entry"]["id"]) == list:
            # remove all non-word characters
            prefix = re.sub(r"[^\w]", " ", self.mmcif["_entry"]["id"][-1])
        else:
            prefix = re.sub(r"[^\w]", " ", self.mmcif["_entry"]["id"])

        # replace all whitespace with a slash
        prefix = re.sub(r"\s+","-", prefix.strip())

        if not in_memory:
            output = "{}.cif".format(prefix)
            print("writing to {}".format(output))
            cifo = CifFileWriter(output)
        else:
            cifo = CifFileWriter(self.output)
        # clean up undefined dictionary
        # force to copy a list to avoid 
        # RuntimeError: dictionary changed size during iteration
        for k in list(self.mmcif): 
            if not self.mmcif[k]:
                del self.mmcif[k]

        cifo.write({ "CCWS" : self.mmcif }) # 'CCWS' ?