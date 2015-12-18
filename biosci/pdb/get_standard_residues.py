from ftplib import FTP
from pprint import pformat

#Log in
ftp = FTP("ftp.wwpdb.org")
ftp.login()
ftp.cwd("pub/pdb/data/monomers")


#Residues required
residues = ["PHE", "TRP", "MET", "ILE", "ASN",
            "THR", "HIS", "GLN", "GLU", "ASP",
            "TYR", "CYS", "ARG", "PRO", "LEU",
            "GLY", "ALA", "VAL", "SER", "LYS"]
residues = {r: {} for r in residues}


#Define line processing function
def process_conect(line):
    chunks = line.split()
    if chunks[0] == "RESIDUE":
        global current_residue
        current_residue = chunks[1]
    elif chunks[0] == "CONECT":
        residues[current_residue][chunks[1]] = chunks[3:]


#Get files
for residue in residues:
    ftp.retrlines("RETR %s" % residue, callback=process_conect)


#Save
f = open("residues.txt", "w")
f.write(pformat(residues))
f.close()
