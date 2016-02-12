from collections import Counter
from .crystal import *
import math
from .exceptions import *

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

        self.models = [Model(d, self.data.miscellaneous.sites,
         self.data.secondary_structure, self.data.connectivity,
          self.data.connectivity_annotation, self.data.heterogen) for d in self.data.coordinates.models]
        self.model = self.models[0]

        self.unit_cell = UnitCell(self.data.crystal)
        self.submission_transformation = SubmittedCoordinatesTransformation(self.data.crystal)
        self.crystal_transformation = CrystallographicCoordinatesTransformation(self.data.crystal)
        self.matrix_transformation = MatrixTransformation(self.data.crystal)



class ChemicalBond:
    """A covalent bond, or similarly strong bond"""

    def __init__(self, *atoms, peptide=False, cis=False, disulphide=False, specified_distance=None):
        assert len(atoms) == 2
        self.atoms = atoms
        self.peptide = peptide
        self.cis = cis
        self.disulphide = disulphide
        self.specified_distance = specified_distance

        atoms[0].bonds.append(self)
        atoms[1].bonds.append(self)


    def __repr__(self):
        return "%s—%s" % (self.atoms[0].element, self.atoms[1].element)



class AtomicStructure:
    """Some structure that contains atoms."""

    def __init__(self, atoms):
        self.atoms = atoms
        if not self.atoms:
            raise PdbStructureError("Structure has no atoms")
        self.mass = sum([a.mass for a in self.atoms])


    def __repr__(self):
        return ", ".join([str(a) for a in self.atoms])


    def get_bonds(self):
        bonds = []
        for atom in self.atoms:
            for bond in atom.bonds:
                if bond not in bonds:
                    bonds.append(bond)
        return bonds


    def get_atom_by_number(self, number):
        for atom in self.atoms:
            if atom.number == number:
                return atom


    def get_atoms_by_name(self, name):
        return [a for a in self.atoms if a.name == name]


    def get_atoms_by_element(self, element):
        return [a for a in self.atoms if a.element == element]


    def get_atom_counts(self):
        atoms = [atom.element.upper() for atom in self.atoms if atom.element.upper() != "H"]
        return Counter(atoms)


    def count_atomic_contacts(self, other_atomic_structure, cutoff):
        """How many atomic interactions are there between this object and
        another atomic structure?"""

        contacts = 0
        for atom in self.atoms:
            nearby = atom.nearby_atoms(cutoff=cutoff, covalent_count=2)
            contacts += len([a for a in nearby if a in other_atomic_structure.atoms and a not in self.atoms])
        return contacts


    def count_internal_atomic_contacts(self, cutoff):
        """How many internal atomic interactions are there in this object?"""

        contacts = 0
        for atom in self.atoms:
            nearby = atom.nearby_atoms(cutoff=cutoff, covalent_count=2)
            contacts += len([a for a in nearby if a in self.atoms])
        return contacts


    def average_x(self):
        """Return the average x coordinate."""

        return sum([atom.x for atom in self.atoms]) / len(self.atoms)


    def average_y(self):
        """Return the average y coordinate."""

        return sum([atom.y for atom in self.atoms]) / len(self.atoms)


    def average_z(self):
        """Return the average z coordinate."""

        return sum([atom.z for atom in self.atoms]) / len(self.atoms)


    def distance_to_other_structure(self, other_atomic_structure):
        """What is the distance between this atomic structure and another, using
        average cartesian coordinates."""

        x_sum = math.pow((other_atomic_structure.average_x() - self.average_x()), 2)
        y_sum = math.pow((other_atomic_structure.average_y() - self.average_y()), 2)
        z_sum = math.pow((other_atomic_structure.average_z() - self.average_z()), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)
        return distance





class ResiduicStructure(AtomicStructure):
    """Some structure that contains residues."""

    def __init__(self, residues):
        self.residues = residues
        atoms = []
        for residue in self.residues:
            atoms += residue.atoms
        AtomicStructure.__init__(self, atoms)


    def __repr__(self):
        return ", ".join([str(r) for r in self.residues])


    def __len__(self):
        return len(self.residues)


    def get_residues_by_chain(self, chain_id):
        return [r for r in self.residues if r.chain.name == chain_id]


    def get_residues_by_name(self, name):
        return [r for r in self.residues if r.name == name]


    def get_residue_by_number(self, number):
        """Returns the first residue with a matching number."""
        for residue in self.residues:
            if residue.number == number:
                return residue


    def get_residues_by_number(self, number):
        """Returns all residues with matching number."""
        return [r for r in self.residues if r.number == number]


    def get_continuous_sequence(self):
        """If all the residues are on the same chain, return an unbroken
        string of residues that contains all residues in this orginal structure."""

        if len(self.residues) > 1:
            if all([r.chain is self.residues[0].chain for r in self.residues[1:]]):
                #All residues are on the same chain
                return ResiduicStructure(
                 [r for r in self.residues[0].chain.residues if
                  r.number >= min([res.number for res in self.residues]) and
                   r.number <= max([res.number for res in self.residues])]
                )
        else:
            return self



class Model(AtomicStructure):
    """A PDB model."""

    def __init__(self, model_dict, site_dicts, secondary_section, connect_section, connect_annotation_section, heterogen_section):
        #Get chains
        chain_ids = sorted(list(set([a["chain_id"] for a in model_dict["atoms"] if not a["het"]])))
        self.chains = [Chain(
         [a for a in model_dict["atoms"] if a["chain_id"] == chain_id and not a["het"]],
         [t for t in model_dict["ters"] if t["chain_id"] == chain_id],
         [h for h in secondary_section.helices if h["start_residue_chain"] == chain_id]
        ) for chain_id in chain_ids]
        atoms = []
        for chain in self.chains:
            atoms += chain.atoms

        #Get hets
        het_ids = sorted(list(set([(a["chain_id"], a["res_seq"]) for a in model_dict["atoms"] if a["het"]])))
        self.hets = [Het(
         [a for a in model_dict["atoms"] if (a["chain_id"], a["res_seq"]) == het_id]
        ) for het_id in het_ids]
        for het in self.hets:
            atoms += het.atoms
            het_dict_matches = [h for h in heterogen_section.hetnams if h["code"] == het.name]
            if het_dict_matches:
                het.full_name = het_dict_matches[0]["fullname"]
            else:
                het.full_name = None
        AtomicStructure.__init__(self, atoms)
        for atom in self.atoms:
            atom.model = self

        #Get sites
        self.pdb_sites = [PdbSite(s, self) for s in site_dicts]
        self.pdb_sites = [site for site in self.pdb_sites if len(site.residues)]

        #Get helices
        self.helices = []
        for chain in self.chains:
            self.helices += chain.helices

        #Get sheets
        self.sheets = [Sheet(s, self) for s in secondary_section.sheets]

        #Connect atoms together (from CONECT records)
        for atom_dict in connect_section.atoms:
            atom_obj = self.get_atom_by_number(atom_dict["atom_id"])
            for bonded_atom_id in atom_dict["bonded_atoms"]:
                bonded_atom_obj = self.get_atom_by_number(bonded_atom_id)
                atom_obj.bond(bonded_atom_obj)

        #Connect atoms together (standard residues)
        from .residues import residues
        for chain in self.chains:
            for residue in chain.residues:
                if residue.name in residues.keys():
                    residue_dict = residues[residue.name]
                    for atom_name in residue_dict.keys():
                        matching_atoms = residue.get_atoms_by_name(atom_name)
                        if len(matching_atoms) == 1:
                            atom = matching_atoms[0]
                            for bonded_atom_name in residue_dict[atom_name]:
                                matching_atoms = residue.get_atoms_by_name(bonded_atom_name)
                                if len(matching_atoms) == 1:
                                    atom.bond(matching_atoms[0])

        #Connect atoms together (peptide bonds)
        for chain in self.chains:
            for index, residue in enumerate(chain.residues[:-1]):
                c, n = None, None
                cs = residue.get_atoms_by_name("C")
                if len(cs) == 1:
                    c = cs[0]
                ns = chain.residues[index+1].get_atoms_by_name("N")
                if len(ns) == 1:
                    n = ns[0]
                if c and n:
                    c.bond(n, peptide=True)

        #Connect atoms together (SS, LINK and CISPEP) (MISSES SOME INFO)
        for ssbond in connect_annotation_section.ssbonds:
            residue1 = self.get_chain_by_name(ssbond["residue_1_chain"]
             ).get_residue_by_number(ssbond["residue_1_number"])
            residue2 = self.get_chain_by_name(ssbond["residue_2_chain"]
             ).get_residue_by_number(ssbond["residue_2_number"])
            atom1 = residue1.get_atoms_by_element("S")[0]
            atom2 = residue2.get_atoms_by_element("S")[0]
            atom1.bond(atom2, disulphide=True, specified_distance=ssbond["disulfide_distance"]) #No symetry information yet
        for link in connect_annotation_section.links:
            atom1 = self.get_atom_by_number(link["residue_1_atom"])
            atom2 = self.get_atom_by_number(link["residue_2_atom"])
            if atom1 and atom2: atom1.bond(atom2)
        for cispep in connect_annotation_section.cispeps:
            residue1 = self.get_chain_by_name(cispep["residue_1_chain"]
             ).get_residue_by_number(cispep["residue_1_number"])
            residue2 = self.get_chain_by_name(cispep["residue_2_chain"]
             ).get_residue_by_number(cispep["residue_2_number"])
            residue1_peptide_bonds = []
            for atom in residue1.atoms:
                for bond in atom.bonds:
                    if bond not in residue1_peptide_bonds and bond.peptide:
                        residue1_peptide_bonds.append(bond)
            residue2_peptide_bonds = []
            for atom in residue2.atoms:
                for bond in atom.bonds:
                    if bond not in residue2_peptide_bonds and bond.peptide:
                        residue2_peptide_bonds.append(bond)
            in_both = [bond for bond in residue1_peptide_bonds if bond in residue2_peptide_bonds]
            if len(in_both) == 1:
                in_both[0].cis = True
                in_both[0].cis_angle = cispep["angle_measure"]




    def __repr__(self):
        return "<Model (%i atoms)>" % len(self.atoms)


    def get_chain_by_name(self, name):
        for chain in self.chains:
            if chain.name == name:
                return chain


    def get_het_by_number(self, number):
        for het in self.hets:
            if het.number == number:
                return het


class Chain(ResiduicStructure):
    "A chain of residues."

    def __init__(self, atoms, termini, helix_dicts):
        self.name = atoms[0]["chain_id"]

        #Get residues
        residue_numbers = sorted(list(set([a["res_seq"] for a in atoms])))
        residues = [Residue(
         [a for a in atoms if a["res_seq"] == residue_number]
        ) for residue_number in residue_numbers]
        ResiduicStructure.__init__(self, residues)
        for residue in self.residues:
            residue.chain = self

        #Specify terminating residue
        if len(termini) == 1:
            matching_residue = [r for r in self.residues if r.number == termini[0]["res_seq"]]
            if matching_residue:
                matching_residue[0].terminus = True

        #Get helices
        self.helices = [Helix(h, self) for h in helix_dicts]

        self.strands = [] #Gets filled in later



    def __repr__(self):
        return "<Chain %s (%i residues)>" % (self.name, len(self.residues))



class Residue(AtomicStructure):
    "An amino acid residue."

    RESIDUE_NAMES = {
     "phenylalanine": ("PHE", "F"), "PHE": ("phenylalanine", "F"), "F": ("phenylalanine", "PHE"),
      "tryptophan": ("TRP", "W"), "TRP": ("tryptophan", "W"), "W": ("tryptophan", "TRP"),
       "methionine": ("MET", "M"), "MET": ("methionine", "M"), "M": ("methionine", "MET"),
        "isoleucine": ("ILE", "I"), "ILE": ("isoleucine", "I"), "I": ("isoleucine", "ILE"),
         "asparagine": ("ASN", "N"), "ASN": ("asparagine", "N"), "N": ("asparagine", "ASN"),
          "threonine": ("THR", "T"), "THR": ("threonine", "T"), "T": ("threonine", "THR"),
           "histidine": ("HIS", "H"), "HIS": ("histidine", "H"), "H": ("histidine", "HIS"),
            "glutamine": ("GLN", "Q"), "GLN": ("glutamine", "Q"), "Q": ("glutamine", "GLN"),
             "glutamate": ("GLU", "E"), "GLU": ("glutamate", "E"), "E": ("glutamate", "GLU"),
              "aspartate": ("ASP", "D"), "ASP": ("aspartate", "D"), "D": ("aspartate", "ASP"),
               "tyrosine": ("TYR", "Y"), "TYR": ("tyrosine", "Y"), "Y": ("tyrosine", "TYR"),
                "cysteine": ("CYS", "C"), "CYS": ("cysteine", "C"), "C": ("cysteine", "CYS"),
                 "arginine": ("ARG", "R"), "ARG": ("arginine", "R"), "R": ("arginine", "ARG"),
                  "proline": ("PRO", "P"), "PRO": ("proline", "P"), "P": ("proline", "PRO"),
                   "leucine": ("LEU", "L"), "LEU": ("leucine", "L"), "L": ("leucine", "LEU"),
                    "glycine": ("GLY", "G"), "GLY": ("glycine", "G"), "G": ("glycine", "GLY"),
                     "alanine": ("ALA", "A"), "ALA": ("alanine", "A"), "A": ("alanine", "ALA"),
                      "valine": ("VAL", "V"), "VAL": ("valine", "V"), "V": ("valine", "VAL"),
                       "serine": ("SER", "S"), "SER": ("serine", "S"), "S": ("serine", "SER"),
                        "lysine": ("LYS", "K"), "LYS": ("lysine", "K"), "K": ("lysine", "LYS")
    }

    def __init__(self, atoms):
        self.number = atoms[0]["res_seq"]
        self.name = atoms[0]["res_name"]
        self.terminus = False

        #Get atoms
        AtomicStructure.__init__(self, [Atom(a) for a in atoms])
        for atom in self.atoms:
            atom.molecule = self


    def __repr__(self):
        return "<%s (%s%i)>" % (self.name, self.chain.name, self.number)


    def connected_residues(self):
        residues = []
        for atom in self.atoms:
            for bonded_atom in atom.bonded_atoms:
                if bonded_atom.molecule is not self and bonded_atom.molecule not in residues:
                    residues.append(bonded_atom.molecule)
        return residues



class Atom:
    """An atom."""

    def __init__(self, atom_dict):
        self.number = atom_dict["serial"]
        self.name = atom_dict["name"]
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
        self.bonds = []


    def __repr__(self):
        return "<%s>" % self.name


    def __getattr__(self, key):
        if key == "bonded_atoms":
            atoms = []
            for bond in self.bonds:
                for atom in bond.atoms:
                    if atom is not self and atom not in atoms:
                        atoms.append(atom)
            return atoms
        else:
            raise AttributeError("Atom object has no attribute %s" % key)


    def bond(self, other_atom, **kwargs):
        if other_atom not in self.bonded_atoms:
            ChemicalBond(self, other_atom, **kwargs)



    def distance_to(self, other_atom):
        x_sum = math.pow((other_atom.x - self.x), 2)
        y_sum = math.pow((other_atom.y - self.y), 2)
        z_sum = math.pow((other_atom.z - self.z), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)
        return distance


    def nearby_atoms(self, cutoff, covalent_count=1):
        atoms_to_exlcude = [self]
        for _ in range(covalent_count):
            for atom in atoms_to_exlcude[:]:
                for bonded_atom in atom.bonded_atoms:
                    if bonded_atom not in atoms_to_exlcude:
                        atoms_to_exlcude.append(bonded_atom)

        atoms = []
        for atom in self.model.atoms:
            if atom not in atoms_to_exlcude and atom.distance_to(self) <= cutoff:
                atoms.append(atom)
        return atoms



class Het(AtomicStructure):
    """A ligand or other non-polymeric molecule (including solvents)."""

    def __init__(self, atoms):
        self.number = atoms[0]["res_seq"]
        self.name = atoms[0]["res_name"]

        #Get atoms
        AtomicStructure.__init__(self, [Atom(a) for a in atoms])
        for atom in self.atoms:
            atom.molecule = self


    def __repr__(self):
        return "<%s (%i atom%s)>" % (self.name, len(self.atoms), "" if len(self.atoms) == 1 else "s")


    def get_nearby_residues(self, cutoff=3):
        """Returns a list of residues close to this ligand."""

        nearby_atoms = []
        for atom in self.atoms:
            for near_atom in atom.nearby_atoms(cutoff):
                if near_atom not in nearby_atoms:
                    nearby_atoms.append(near_atom)
        nearby_atoms = [atom for atom in nearby_atoms if atom not in self.atoms
         and isinstance(atom.molecule, Residue)]
        residues = list(set([atom.molecule for atom in nearby_atoms]))
        residues = sorted(residues, key=lambda k: k.chain.name)
        residues = sorted(residues, key=lambda k: k.number)
        return residues



class PdbSite(ResiduicStructure):
    """A site named in the PDB file."""

    def __init__(self, site_dict, model):
        self.name = site_dict["name"]

        residues = []
        self.hets = []
        for residue_dict in site_dict["residues"]:
            chain = model.get_chain_by_name(residue_dict["chain"])
            if chain:
                residue = chain.get_residue_by_number(residue_dict["residue_number"])
                if residue:
                    residues.append(residue)
            het = model.get_het_by_number(residue_dict["residue_number"])
            if het:
                self.hets.append(het)
        ResiduicStructure.__init__(self, residues)


    def __repr__(self):
        return "<%s: %s>" % (self.name, ", ".join([str(r) for r in self.residues + self.hets]))


class Helix(ResiduicStructure):
    """An alpha helix."""

    CLASSES = {1: "Right-handed alpha", 2: "Right-handed omega", 3: "Right-handed pi",
     4: "Right-handed gamma", 5: "Right-handed 3 - 10", 6: "Left-handed alpha", 7: "Left-handed omega",
      8: "Left-handed gamma", 9: "2 - 7 ribbon/helix", 10: "Polyproline"}

    def __init__(self, helix_dict, chain):
        self.number = helix_dict["serial"]
        self.name = helix_dict["helix_id"]
        self.chain = chain
        residues = [chain.get_residue_by_number(i) for i in
         range(helix_dict["start_residue_number"], helix_dict["end_residue_number"] + 1)]
        residues = [r for r in residues if r] #Important if file lists residue that doesn't exist
        ResiduicStructure.__init__(self, residues)
        self.helix_class = self.CLASSES[helix_dict["helix_class"] if helix_dict["helix_class"] else 1]
        self.comment = helix_dict["comment"]



class Sheet(ResiduicStructure):
    """A beta sheet."""

    def __init__(self, sheet_dict, model):
        self.name = sheet_dict["sheet_id"]
        self.strands = [Strand(sheet_dict["strands"][0], model)]
        for strand_dict in sheet_dict["strands"][1:]:
            self.strands.append(Strand(strand_dict, model, self.strands[-1]))
        residues = []
        for strand in self.strands:
            residues += strand.residues
        ResiduicStructure.__init__(self, residues)


    def __repr__(self):
        return "<Beta sheet %s (%i strands)>" % (self.name, len(self.strands))



class Strand(ResiduicStructure):
    """A beta strand."""

    def __init__(self, strand_dict, model, previous_strand=None):
        self.number = strand_dict["strand_id"]
        self.name = self.number
        self.chain = model.get_chain_by_name(strand_dict["start_residue_chain"])
        self.chain.strands.append(self)
        residues = [self.chain.get_residue_by_number(i) for i in
         range(strand_dict["start_residue_number"], strand_dict["end_residue_number"] + 1)]
        residues = [r for r in residues if r] #Important if file lists residue that doesn't exist
        ResiduicStructure.__init__(self, residues)
        self.sense = strand_dict["sense"]

        if previous_strand:
            self.previous_strand = previous_strand
            self.registration = (
             self.get_atoms_by_name(strand_dict["reg_cur_atom"])[0]
              if self.get_atoms_by_name(strand_dict["reg_cur_atom"]) else None,
             self.get_atoms_by_name(strand_dict["reg_prev_atom"])[0]
              if self.get_atoms_by_name(strand_dict["reg_prev_atom"]) else None
            )
        else:
            self.previous_strand, self.registration = None, None
