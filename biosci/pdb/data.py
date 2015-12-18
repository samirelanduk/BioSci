from .exceptions import *

class PdbDataStructure:
    """A processed PdbFile, mostly in dictionary form"""

    def __init__(self, pdb_file):
        self.file = pdb_file

        #Create the sections
        self.secondary_structure = SecondaryStructureSection(self.file)
        self.miscellaneous = MiscellaneousSection(self.file)
        self.crystal = CrystalSection(self.file)
        self.coordinates = CoordinateSection(self.file)
        self.connectivity = ConnectivitySection(self.file)
        self.bookkeeping = BookkeepingSection(self.file)



class PdbSection:

    RECORD_NAMES = ()

    def __init__(self, pdb_file):
        self.records = [r for r in pdb_file.records if r.name in self.RECORD_NAMES]


    def get_records_by_name(self, name):
        return [r for r in self.records if r.name.upper() == name.upper()]



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
                          line[(11 * (x+1)) + 1: (11 * (x+1)) + 5].strip() else None,
                         "i_code": line[((x+1) * 11) + 5] if line[((x+1) * 11) + 5].strip() else None,
                        })
            self.sites.append(site)



class CrystalSection(PdbSection):

    RECORD_NAMES = ("CRYST1", "MTRIX1", "MTRIX2", "MTRIX3",
     "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process CRYST1
        cryst1 = [r for r in self.records if r.name == "CRYST1"]
        if len(cryst1) == 0:
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.s_group, self.z = (
             None, None, None, None, None, None, None, None)
        else:
            cryst1 = cryst1[0]
            self.a = float(cryst1[6:15].strip()) if cryst1[6:15].strip() else None
            self.b = float(cryst1[15:24].strip()) if cryst1[15:24].strip() else None
            self.c = float(cryst1[24:33].strip()) if cryst1[24:33].strip() else None
            self.alpha = float(cryst1[33:40].strip()) if cryst1[33:40].strip() else None
            self.beta = float(cryst1[40:47].strip()) if cryst1[40:47].strip() else None
            self.gamma = float(cryst1[47:54].strip()) if cryst1[47:54].strip() else None
            self.s_group = cryst1[55:66].strip() if cryst1[55:66].strip() else None
            self.z = int(cryst1[66:70].strip()) if cryst1[66:70].strip() else None

        #Process ORIGXns
        origx1 = [r for r in self.records if r.name == "ORIGX1"]
        origx2 = [r for r in self.records if r.name == "ORIGX2"]
        origx3 = [r for r in self.records if r.name == "ORIGX3"]
        if len(origx1) == 0:
            self.o11, self.o12, self.o13, self.t1 = None, None, None, None
        else:
            origx1 = origx1[0]
            self.o11 = float(origx1[10:20].strip()) if origx1[10:20].strip() else None
            self.o12 = float(origx1[20:30].strip()) if origx1[20:30].strip() else None
            self.o13 = float(origx1[30:40].strip()) if origx1[30:40].strip() else None
            self.t1 = float(origx1[45:55].strip()) if origx1[45:55].strip() else None
        if len(origx2) == 0:
            self.o21, self.o22, self.o23, self.t2 = None, None, None, None
        else:
            origx2 = origx2[0]
            self.o21 = float(origx2[10:20].strip()) if origx2[10:20].strip() else None
            self.o22 = float(origx2[20:30].strip()) if origx2[20:30].strip() else None
            self.o23 = float(origx2[30:40].strip()) if origx2[30:40].strip() else None
            self.t2 = float(origx2[45:55].strip()) if origx2[45:55].strip() else None
        if len(origx3) == 0:
            self.o31, self.o32, self.o33, self.t3 = None, None, None, None
        else:
            origx3 = origx3[0]
            self.o31 = float(origx3[10:20].strip()) if origx3[10:20].strip() else None
            self.o32 = float(origx3[20:30].strip()) if origx3[20:30].strip() else None
            self.o33 = float(origx3[30:40].strip()) if origx3[30:40].strip() else None
            self.t3 = float(origx3[45:55].strip()) if origx3[45:55].strip() else None

        #Process SCALEns
        scale1 = [r for r in self.records if r.name == "SCALE1"]
        scale2 = [r for r in self.records if r.name == "SCALE2"]
        scale3 = [r for r in self.records if r.name == "SCALE3"]
        if len(scale1) == 0:
            self.s11, self.s12, self.s13, self.u1 = None, None, None, None
        else:
            scale1 = scale1[0]
            self.s11 = float(scale1[10:20].strip()) if scale1[10:20].strip() else None
            self.s12 = float(scale1[20:30].strip()) if scale1[20:30].strip() else None
            self.s13 = float(scale1[30:40].strip()) if scale1[30:40].strip() else None
            self.u1 = float(scale1[45:55].strip()) if scale1[45:55].strip() else None
        if len(scale2) == 0:
            self.s21, self.s22, self.s23, self.u2 = None, None, None, None
        else:
            scale2 = scale2[0]
            self.s21 = float(scale2[10:20].strip()) if scale2[10:20].strip() else None
            self.s22 = float(scale2[20:30].strip()) if scale2[20:30].strip() else None
            self.s23 = float(scale2[30:40].strip()) if scale2[30:40].strip() else None
            self.u2 = float(scale2[45:55].strip()) if scale2[45:55].strip() else None
        if len(scale3) == 0:
            self.s31, self.s32, self.s33, self.u3 = None, None, None, None
        else:
            scale3 = scale3[0]
            self.s31 = float(scale3[10:20].strip()) if scale3[10:20].strip() else None
            self.s32 = float(scale3[20:30].strip()) if scale3[20:30].strip() else None
            self.s33 = float(scale3[30:40].strip()) if scale3[30:40].strip() else None
            self.u3 = float(scale3[45:55].strip()) if scale3[45:55].strip() else None

        #Process MTRIXns
        mtrix1 = [r for r in self.records if r.name == "MTRIX1"]
        mtrix2 = [r for r in self.records if r.name == "MTRIX2"]
        mtrix3 = [r for r in self.records if r.name == "MTRIX3"]
        if len(mtrix1) == 0:
            self.m11, self.m12, self.m13, self.v1 = None, None, None, None
        else:
            mtrix1 = mtrix1[0]
            self.m11 = float(mtrix1[10:20].strip()) if mtrix1[10:20].strip() else None
            self.m12 = float(mtrix1[20:30].strip()) if mtrix1[20:30].strip() else None
            self.m13 = float(mtrix1[30:40].strip()) if mtrix1[30:40].strip() else None
            self.v1 = float(mtrix1[45:55].strip()) if mtrix1[45:55].strip() else None
        if len(mtrix2) == 0:
            self.m21, self.m22, self.m23, self.v2 = None, None, None, None
        else:
            mtrix2 = mtrix2[0]
            self.m21 = float(mtrix2[10:20].strip()) if mtrix2[10:20].strip() else None
            self.m22 = float(mtrix2[20:30].strip()) if mtrix2[20:30].strip() else None
            self.m23 = float(mtrix2[30:40].strip()) if mtrix2[30:40].strip() else None
            self.v2 = float(mtrix2[45:55].strip()) if mtrix2[45:55].strip() else None
        if len(mtrix3) == 0:
            self.m31, self.m32, self.m33, self.v3 = None, None, None, None
        else:
            mtrix3 = mtrix3[0]
            self.m31 = float(mtrix3[10:20].strip()) if mtrix3[10:20].strip() else None
            self.m32 = float(mtrix3[20:30].strip()) if mtrix3[20:30].strip() else None
            self.m33 = float(mtrix3[30:40].strip()) if mtrix3[30:40].strip() else None
            self.v3 = float(mtrix3[45:55].strip()) if mtrix3[45:55].strip() else None



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


class ConnectivitySection(PdbSection):

    RECORD_NAMES = ("CONECT")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process CONECTs
        atoms = list(set([int(r[6:11].strip()) for r in self.records]))
        self.atoms = []
        for atom_id in atoms:
            atom = {"atom_id": atom_id, "bonded_atoms": []}
            strings = [r[11:31] for r in self.records if r[6:11].strip() == str(atom_id)]
            for string in strings:
                if string[:5].strip():
                    atom["bonded_atoms"].append(int(string[:5].strip()))
                if string[5:10].strip():
                    atom["bonded_atoms"].append(int(string[5:10].strip()))
                if string[10:15].strip():
                    atom["bonded_atoms"].append(int(string[10:15].strip()))
                if string[15:20].strip():
                    atom["bonded_atoms"].append(int(string[15:20].strip()))
            self.atoms.append(atom)



class BookkeepingSection(PdbSection):

    RECORD_NAMES = ("MASTER", "END")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process MASTER
        master = self.get_records_by_name("MASTER")
        if master:
            master = master[0]
            self.num_remark = int(master[10:15].strip()) if master[10:15].strip() else None
            self.num_het = int(master[20:25].strip()) if master[20:25].strip() else None
            self.num_helix = int(master[25:30].strip()) if master[25:30].strip() else None
            self.num_sheet = int(master[30:35].strip()) if master[30:35].strip() else None
            self.num_site = int(master[40:45].strip()) if master[40:45].strip() else None
            self.num_xform = int(master[45:50].strip()) if master[45:50].strip() else None
            self.num_coord = int(master[50:55].strip()) if master[50:55].strip() else None
            self.num_ter = int(master[55:60].strip()) if master[55:60].strip() else None
            self.num_conect = int(master[60:65].strip()) if master[60:65].strip() else None
            self.num_seq = int(master[65:70].strip()) if master[65:70].strip() else None
        else:
            self.num_remark = None
            self.num_het = None
            self.num_helix = None
            self.num_sheet = None
            self.num_site = None
            self.num_xform = None
            self.num_coord = None
            self.num_ter = None
            self.num_conect = None
            self.num_seq = None

        #Process END
        if not self.get_records_by_name("END"):
            raise PdbDataError("This PDB has no END")
