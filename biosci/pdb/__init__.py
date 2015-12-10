from .file import *
from .data import *
from .structure import *
from .exceptions import *
import requests

def get_from_code(code):
    """Takes a 4-char PDB identifier and gets the PDB over the internet,
     then processes it."""

    response = requests.get(
     "http://www.rcsb.org/pdb/files/%s.pdb" % code
    )
    if response.status_code == 200 and response.text[:6] == "HEADER":
        contents = response.text
    else:
        raise PdbError("%s does not seem to be a valid PDB code." % code)

    pdb_file = PdbFile(contents)
    pdb_data = PdbDataStructure(pdb_file)
    return PdbStructure(pdb_data)
