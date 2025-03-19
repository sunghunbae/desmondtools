from pathlib import Path
from dotmap import DotMap
from desmondtools import Multisim



class Event:
    def __init__(self, filename:Path | str):
        if isinstance(filename, Path):
            result = Multisim.expr.parse_file(filename.as_posix())
        else:
            result = Multisim.expr.parse_file(filename)
        d = result.as_dict()
        Multisim.traverse_dict(d)
        self.dot = DotMap(d)
        


class Interaction(Event):
    def __init__(self, filename:Path | str, verbose:bool=False):
        super().__init__(filename)

        self.num_frames = 0
        self.HBond = {}
        self.Hydrophobic = {}
        self.WaterBridge = {}
        self.Polar = {}
        self.HalogenBond = {}
        self.LigWat = {}
        self.Metal = {}
        self.PiCat = {}
        self.PiPi = {}

        for section in self.dot.Keywords:
            try:
                assert section.ProtLigInter.HBondResult
                self.num_frames = len(section.ProtLigInter.HBondResult)
                for frame in section.ProtLigInter.HBondResult:
                    # [[3 "_:ARG_143:HH22" d-s "L-FRAG_0:N6" ]]
                    for (frameno, prot, hbond_type, lig) in frame:
                        prot = prot.strip('\"')
                        (_, resid, atom) = prot.split(":")
                        (resName, resSeq) = resid.split("_")
                        resSeq = int(resSeq)
                        if resSeq in self.HBond:
                            self.HBond[resSeq]['count'] += 1
                        else:
                            self.HBond[resSeq] = {'resName': resName, 'count':1 }
                for resSeq in sorted(self.HBond):
                    fraction = float(self.HBond[resSeq]['count'])/self.num_frames
                    if verbose:
                        print(f"HBond {self.HBond[resSeq]['resName']}_{resSeq} {fraction:5.3f} {self.num_frames}")
            except:
                pass

            try:
                assert section.ProtLigInter.HydrophobicResult
                self.num_frames = len(section.ProtLigInter.HydrophobicResult)
                for frame in section.ProtLigInter.HydrophobicResult:
                    # [[0 "_:PHE_223" L-FRAG_0 ] [0 "_:ALA_241" L-FRAG_0 ]]
                    for (frameno, prot, lig) in frame:
                        prot = prot.strip('\"')
                        (_, resid) = prot.split(":")
                        (resName, resSeq) = resid.split("_")
                        resSeq = int(resSeq)
                        if resSeq in self.Hydrophobic:
                            self.Hydrophobic[resSeq]['count'] += 1
                        else:
                            self.Hydrophobic[resSeq] = {'resName': resName, 'count':1 }
                for resSeq in sorted(self.Hydrophobic):
                    fraction = float(self.Hydrophobic[resSeq]['count'])/self.num_frames
                    if verbose:
                        print(f"Hydrophobic {self.Hydrophobic[resSeq]['resName']}_{resSeq} {fraction:5.3f} {self.num_frames}")
            except:
                pass
            
            try:
                assert section.ProtLigInter.PolarResult
                self.num_frames = len(section.ProtLigInter.PolarResult)
                for frame in section.ProtLigInter.PolarResult:
                    # [[1 "_:GLU_216:OE2" b "L-FRAG_1:N3" 4.45 ]]
                    for (frameno, prot, _, lig, _) in frame:
                        prot = prot.strip('\"')
                        (_, resid, atom) = prot.split(":")
                        (resName, resSeq) = resid.split("_")
                        resSeq = int(resSeq)
                        if resSeq in self.Polar:
                            self.Polar[resSeq]['count'] += 1
                        else:
                            self.Polar[resSeq] = {'resName': resName, 'count':1 }
                for resSeq in sorted(self.Polar):
                    fraction = float(self.Polar[resSeq]['count'])/self.num_frames
                    if verbose:
                        print(f"Polar {self.Polar[resSeq]['resName']}_{resSeq} {fraction:5.3f} {self.num_frames}")
            except:
                pass

            try:
                assert section.ProtLigInter.WaterBridgeResult
                self.num_frames = len(section.ProtLigInter.WaterBridgeResult)
                for frame in section.ProtLigInter.WaterBridgeResult:
                    # [[3 "_:GLU_216:OE2" a "L-FRAG_0:N2" a 2431 ]]
                    for (frameno, prot, _, lig, _, _) in frame:
                        prot = prot.strip('\"')
                        (_, resid, atom) = prot.split(":")
                        (resName, resSeq) = resid.split("_")
                        resSeq = int(resSeq)
                        if resSeq in self.WaterBridge:
                            self.WaterBridge[resSeq]['count'] += 1
                        else:
                            self.WaterBridge[resSeq] = {'resName': resName, 'count':1 }
                for resSeq in sorted(self.WaterBridge):
                    fraction = float(self.WaterBridge[resSeq]['count'])/self.num_frames
                    if verbose:
                        print(f"WaterBridge {self.WaterBridge[resSeq]['resName']}_{resSeq} {fraction:5.3f} {self.num_frames}")
            except:
                pass