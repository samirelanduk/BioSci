from .exceptions import *

class PdbDataStructure:
    """A processed PdbFile, mostly in dictionary form"""

    def __init__(self, pdb_file):
        self.file = pdb_file

        #Create the sections
        self.coordinates = CoordinateSection(self.file)
        self.secondary_structure = SecondaryStructureSection(self.file)
        self.miscellaneous = MiscellaneousSection(self.file)



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
            models = []
            for rec in model_recs:
                end = [e for e in self.get_records_by_name("ENDMDL") if e.number > rec.number][0]
                model = [r for r in self.records if r.number > rec.number and r.number < end.number]
                models.append(model)
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
                matching_atom = [a for a in model["atoms"] if a["serial"] == int(anisou[6:11].strip())]
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



class SecondaryStructureSection(PdbSection):

    RECORD_NAMES = ("HELIX", "SHEET")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process HELIXs
        helices = [r for r in self.records if r.name == "HELIX"]
        self.helices = [{
         "serial": int(l[7:10].strip()) if l[7:10].strip() else None,
         "helix_id": l[11:14].strip() if l[11:14].strip() else None,
         "start_residue_name": l[15:18].strip() if l[15:18].strip() else None,
         "start_residue_chain": l[19] if l[19].strip() else None,
         "start_residue_number": int(l[21:25].strip()) if l[21:25].strip() else None,
         "start_residue_insert": l[25] if l[25].strip() else None,
         "end_residue_name": l[27:30].strip() if l[27:28].strip() else None,
         "end_residue_chain": l[31] if l[31].strip() else None,
         "end_residue_number": int(l[33:37].strip()) if l[33:37].strip() else None,
         "end_residue_insert": l[37] if l[37].strip() else None,
         "helix_class": int(l[38:40].strip()) if l[38:40].strip() else None,
         "comment": l[40:70].strip() if l[40:70].strip() else None,
         "length": int(l[71:76].strip()) if l[71:76].strip() else None
        } for l in helices]

        #Process SHEETs
        sheet_lines = [r for r in self.records if r.name == "SHEET"]
        sheets = list(set([l[11:14].strip() for l in sheet_lines]))
        self.sheets = []
        for sheet in sheets:
            lines = [l for l in sheet_lines if l[11:14].strip() == sheet]
            strands = [{
             "strand_id": int(l[7:10].strip()) if l[7:10].strip() else None,
             "start_residue_name": l[17:20].strip() if l[17:20].strip() else None,
             "start_residue_chain": l[21] if l[21].strip() else None,
             "start_residue_number": int(l[22:26].strip()) if l[22:26].strip() else None,
             "start_residue_insert": l[26] if l[26].strip() else None,
             "end_residue_name": l[28:31].strip() if l[28:31].strip() else None,
             "end_residue_chain": l[32] if l[32].strip() else None,
             "end_residue_number": int(l[33:37].strip()) if l[33:37].strip() else None,
             "end_residue_insert": l[37] if l[37].strip() else None,
             "sense": int(l[38:40].strip()) if l[38:40].strip() else None,
             "reg_cur_atom": l[41:45].strip() if l[41:45].strip() else None,
             "reg_cur_residue": l[45:48] if l[45:48].strip() else None,
             "reg_cur_chain": l[49] if l[49].strip() else None,
             "reg_cur_number": int(l[50:54].strip()) if l[50:54].strip() else None,
             "reg_cur_insert": l[54] if l[54].strip() else None,
             "reg_prev_atom": l[56:60].strip() if l[56:60].strip() else None,
             "reg_prev_residue": l[60:63] if l[60:63].strip() else None,
             "reg_prev_chain": l[64] if l[64].strip() else None,
             "reg_prev_number": int(l[65:69].strip()) if l[65:69].strip() else None,
             "reg_prev_insert": l[69] if l[69] else None
            } for l in lines]
            self.sheets.append({"sheet_id":sheet, "strands":strands})



class MiscellaneousSection(PdbSection):

    RECORD_NAMES = ("SITE")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process SITEs
        sites = [r for r in self.records if r.name == "SITE"]
        site_names = sorted(list(set([s[11:14].strip() for s in sites])))
        self.sites = []
        for name in site_names:
            lines = [r for r in sites if r[11:14].strip() == name]
            site = {"name": name, "residues": []}
            for line in lines:
                for x in (1,2,3,4):
                    if line[(11 * x) + 7: (11 * x) + 16].strip():
                        site["residues"].append({
                         "chain": line[(x+1) * 11] if line[(x+1) * 11].strip() else None,
                         "residue_name": line[(11 * x) + 7: (11 * x) + 10].strip()
                          if line[(11 * x) + 7: (11 * x) + 10].strip() else None,
                         "residue_number": int(line[(11 * (x+1)) + 1: (11 * (x+1)) + 5].strip()) if
                          line[(11 * (x+1)) + 1: (11 * (x+1)) + 5].strip() else None
                        })
            self.sites.append(site)
