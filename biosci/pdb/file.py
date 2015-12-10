from .exceptions import *

VALID_RECORDS = ("HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND", "SOURCE", "KEYWDS",
 "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK",
  "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES",
   "HET", "FORMUL", "HETNAM", "HETSYN",
    "HELIX", "SHEET",
     "SSBOND", "LINK", "CISPEP",
      "SITE",
       "CRYST1", "MTRIX1", "MTRIX2", "MTRIX3", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3",
        "MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL",
         "CONECT",
          "MASTER", "END",
           "HYDBND", "SLTBRG")

class PdbFile:
    """A representation of the PDB file itself, not the structure it represents."""

    def __init__(self, pdb_contents):
        lines = [line for line in pdb_contents.split("\n") if line.strip()]
        self.records = [Record(index, line) for index, line in enumerate(lines, start=1)]


    def get_records_by_name(self, name):
        return [r for r in self.records if r.name.upper() == name.upper()]


    def to_file_contents(self, delim="\n"):
        """Get text for saving to file."""
        return delim.join([r.text for r in self.records])



class Record:
    """A PDB record (a line in the file)."""

    def __init__(self, number, text):
        self.number = number
        self.text = text + (" " * (80 - len(text)))

        #What is the name, and is it valid?
        self.name = text[:6].strip()
        if self.name not in VALID_RECORDS:
            raise PdbFileError("%s is not a valid Record name" % self.name)


    def __repr__(self):
        return "<%s record>" % self.name


    def __getitem__(self, key):
        return self.text[key]
