PERIODIC_TABLE = {
 "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811, "C": 12.0107,
  "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797, "NA": 22.9897, "MG": 24.305,
   "AL": 26.9815, "SI": 28.0855, "P": 30.9738, "S": 32.065, "CL": 35.453, "K": 39.0983,
    "AR": 39.948, "CA": 40.078, "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961,
     "MN": 54.938, "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
      "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904, "KR": 83.8,
       "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224, "NB": 92.9064, "MO": 95.94,
        "TC": 98, "RU": 101.07, "RH": 102.9055, "PD": 106.42, "AG": 107.8682, "CD": 112.411,
         "IN": 114.818, "SN": 118.71, "SB": 121.76, "I": 126.9045, "TE": 127.6,
          "XE": 131.293, "CS": 132.9055, "BA": 137.327, "LA": 138.9055, "CE": 140.116,
           "PR": 140.9077, "ND": 144.24, "PM": 145, "SM": 150.36, "EU": 151.964,
            "GD": 157.25, "TB": 158.9253, "DY": 162.5, "HO": 164.9303, "ER": 167.259,
             "TM": 168.9342, "YB": 173.04, "LU": 174.967, "HF": 178.49, "TA": 180.9479,
              "W": 183.84, "RE": 186.207, "OS": 190.23, "IR": 192.217, "PT": 195.078,
               "AU": 196.9665, "HG": 200.59, "TL": 204.3833, "PB": 207.2, "BI": 208.9804,
                "PO": 209, "AT": 210, "RN": 222, "FR": 223, "RA": 226, "AC": 227,
                 "PA": 231.0359, "TH": 232.0381, "NP": 237, "U": 238.0289, "AM": 243,
                  "PU": 244, "CM": 247, "BK": 247, "CF": 251, "ES": 252, "FM": 257,
                   "MD": 258, "NO": 259, "RF": 261, "LR": 262, "DB": 262, "BH": 264,
                    "SG": 266, "MT": 268, "RG": 272, "HS": 277, "X": 0, "D": 0}

class PdbStructure:
    """A representation of the contents of a PDB file."""

    def __init__(self, pdb_data):
        self.data = pdb_data

        self.models = [Model(d) for d in self.data.coordinates.models]
        self.model = self.models[0]



class Model:
    """A PDB model."""

    def __init__(self, model_dict):
        #Get chains
        chain_ids = sorted(list(set([a["chain_id"] for a in model_dict["atoms"] if not a["het"]])))
        self.chains = [Chain(
         [a for a in model_dict["atoms"] if a["chain_id"] == chain_id and not a["het"]],
         [t for t in model_dict["ters"] if t["chain_id"] == chain_id]
        ) for chain_id in chain_ids]
        self.mass = sum([c.mass for c in self.chains])

        #Get hets
        het_numbers = sorted(list(set([a["res_seq"] for a in model_dict["atoms"] if a["het"]])))
        self.hets = [Het(
         [a for a in model_dict["atoms"] if a["res_seq"] == het_number]
        ) for het_number in het_numbers]
        self.mass += sum([h.mass for h in self.hets])



class Chain:
    "A chain of residues."

    def __init__(self, atoms, termini):
        self.name = atoms[0]["chain_id"]

        #Get residues
        residue_numbers = sorted(list(set([a["res_seq"] for a in atoms])))
        self.residues = [Residue(
         [a for a in atoms if a["res_seq"] == residue_number]
        ) for residue_number in residue_numbers]
        for residue in self.residues:
            residue.chain = self
        self.mass = sum([r.mass for r in self.residues])

        #Specify terminating residue
        if len(termini) == 1:
            matching_residue = [r for r in self.residues if r.number == termini[0]["res_seq"]]
            if matching_residue:
                matching_residue[0].terminus = True


    def __repr__(self):
        return "<Chain %s>" % self.name



class Residue:
    "An amino acid residue."

    def __init__(self, atoms):
        self.number = atoms[0]["res_seq"]
        self.name = atoms[0]["res_name"]
        self.terminus = False

        #Get atoms
        self.atoms = [Atom(a) for a in atoms]
        for atom in self.atoms:
            atom.molecule = self
        self.mass = sum([a.mass for a in self.atoms])


    def __repr__(self):
        return "<%s (%i)>" % (self.name, self.number)



class Atom:
    """An atom."""

    def __init__(self, atom_dict):
        self.number = atom_dict["serial"]
        self.name = atom_dict["serial"]
        self.insert_code = atom_dict["i_code"]
        self.x = atom_dict["x"]
        self.y = atom_dict["y"]
        self.z = atom_dict["z"]
        self.occupancy = atom_dict["occupancy"]
        self.temp_factor = atom_dict["temp_factor"]
        self.element = atom_dict["element"]
        self.charge = atom_dict["charge"]
        self.u11 = atom_dict.get("u11", None)
        self.u22 = atom_dict.get("u22", None)
        self.u33 = atom_dict.get("u33", None)
        self.u12 = atom_dict.get("u12", None)
        self.u13 = atom_dict.get("u13", None)
        self.u23 = atom_dict.get("u23", None)

        self.mass = PERIODIC_TABLE[self.element.upper()]



class Het(Residue):

    def __repr__(self):
        return "<%s (%i atom%s)>" % (self.name, len(self.atoms), "" if len(self.atoms) == 1 else "s")
