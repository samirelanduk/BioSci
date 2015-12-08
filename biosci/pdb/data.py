from .exceptions import *

class PdbDataStructure:
    """A processed PdbFile, mostly in dictionary form"""

    def __init__(self, pdb_file):
        self.file = pdb_file

        #Create the sections
        self.coordinates = CoordinateSection(self.file)



class PdbSection:

    RECORD_NAMES = ()

    def __init__(self, pdb_file):
        self.records = [r for r in pdb_file.records if r.name in self.RECORD_NAMES]


    def get_records_by_name(self, name):
        return [r for r in self.records if r.name.upper() == name.upper()]



class CoordinateSection(PdbSection):

    RECORD_NAMES = ("MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process MODELs
        model_recs = self.get_records_by_name("MODEL")
        if model_recs:
            #There are seperate models specified in the PDB file
            if len(model_recs) != len(self.get_records_by_name("ENDMDL")):
                raise PdbDataError
            models = [[r for r in self.records if
             r.number > model_rec.number and r.number <
              [e for e in self.get_records_by_name("ENDMDL") if e.number > r.number][0].number]
               for model_rec in model_recs]
        else:
            #There is only one one model
            models = [self.records]
        self.models = []

        for model_lines in models:
            model = {}

            #Process ATOMs
            atoms = [r for r in model_lines if r.name == "ATOM"]
            model["atoms"] = [{
             "het": False,
             "serial": int(a[6:11].strip()) if a[6:11].strip() else None,
             "name": a[12:16].strip() if a[12:16].strip() else None,
             "alt_loc": a[16] if a[16].strip() else None,
             "res_name": a[17:20].strip() if a[17:20].strip() else None,
             "chain_id": a[21] if a[21].strip() else None,
             "res_seq": int(a[22:26].strip()) if a[22:26].strip() else None,
             "i_code": a[26] if a[26].strip() else None,
             "x": float(a[30:38].strip()) if a[30:38].strip() else None,
             "y": float(a[38:46].strip()) if a[38:46].strip() else None,
             "z": float(a[46:54].strip()) if a[46:54].strip() else None,
             "occupancy": float(a[54:60].strip()) if a[54:60].strip() else None,
             "temp_factor": float(a[60:66].strip()) if a[60:66].strip() else None,
             "element": a[76:78].strip() if a[76:78].strip() else None,
             "charge": a[78:80].strip() if a[78:80].strip() else None
            } for a in atoms]

            #Process HETATMs
            het_atoms = [r for r in model_lines if r.name == "HETATM"]
            model["atoms"] += [{
             "het": True,
             "serial": int(a[6:11].strip()) if a[6:11].strip() else None,
             "name": a[12:16].strip() if a[12:16].strip() else None,
             "alt_loc": a[16] if a[16].strip() else None,
             "res_name": a[17:20].strip() if a[17:20].strip() else None,
             "chain_id": a[21] if a[21].strip() else None,
             "res_seq": int(a[22:26].strip()) if a[22:26].strip() else None,
             "i_code": a[26] if a[26].strip() else None,
             "x": float(a[30:38].strip()) if a[30:38].strip() else None,
             "y": float(a[38:46].strip()) if a[38:46].strip() else None,
             "z": float(a[46:54].strip()) if a[46:54].strip() else None,
             "occupancy": float(a[54:60].strip()) if a[54:60].strip() else None,
             "temp_factor": float(a[60:66].strip()) if a[60:66].strip() else None,
             "element": a[76:78].strip() if a[76:78].strip() else None,
             "charge": a[78:80].strip() if a[78:80].strip() else None
            } for a in het_atoms]

            #Process ANISOUs
            anisous = [r for r in model_lines if r.name == "ANISOU"]
            for anisou in anisous:
                matching_atom = [a for a in model["atoms"] if atom["serial"] == int(a[6:11].strip())]
                if matching_atom:
                    matching_atom[0]["u11"] = int(anisou[28:35].strip()) if anisou[28:35].strip() else None
                    matching_atom[0]["u22"] = int(anisou[35:42].strip()) if anisou[35:42].strip() else None
                    matching_atom[0]["u33"] = int(anisou[42:49].strip()) if anisou[42:49].strip() else None
                    matching_atom[0]["u12"] = int(anisou[49:56].strip()) if anisou[49:56].strip() else None
                    matching_atom[0]["u13"] = int(anisou[56:63].strip()) if anisou[56:63].strip() else None
                    matching_atom[0]["u23"] = int(anisou[63:70].strip()) if anisou[63:70].strip() else None

            #Process TERs
            ters = [r for r in model_lines if r.name == "TER"]
            model["ters"] = [{
             "serial": int(t[6:11].strip()) if t[6:11].strip() else None,
             "res_name": t[17:20].strip() if t[17:20].strip() else None,
             "chain_id": t[21] if t[21].strip() else None,
             "res_seq": int(t[22:26].strip()) if t[22:26].strip() else None,
             "i_code": t[26] if t[26].strip() else None,
            } for t in ters]

            self.models.append(model)
