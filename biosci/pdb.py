from __future__ import print_function
import re
import math

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


def break_into_sections(sequence, break_points):
    #Always start from the beginning
    if break_points[0] != 0:
        break_points = [0] + break_points

    #Get sections
    sections = []
    for index in break_points:
        if index == break_points[-1]:
            sections.append(sequence[index:])
        else:
            sections.append(sequence[index:break_points[break_points.index(index)+1]])

    return sections



def get_substring_positions(string, substring):
    index = 0
    positions = []

    while substring in string[index:]:
        positions.append(string[index:].find(substring) + index)
        index = positions[-1] + 1

    return positions



class NotValidPdbError(Exception):
    pass




class Record:

    def __init__(self, line, number):
        self.text = line
        self.line_num = number

        #What is the name, and is it valid?
        self.name = line[:6].strip()
        valid = False
        for rec in Pdb.VALID_RECORDS:
            if re.match(rec, self.name):
                valid=True
                break
        if not valid:
            raise NotValidPdbError("%s is not a valid Record (line %i)" %
             (self.name, self.line_num))



    def __repr__(self):
        return self.text.strip()



    def __getitem__(self, key):
        return self.text[key]




class Section:

    def __init__(self, records):
        self.records = records



    def __contains__(self, key):
        return key in self.records



    def __repr__(self):
        return "\n".join([str(r) for r in self.records])


    def _get_record_contents(self, record_name, delim=" ", line_start=10):
        lines = [r for r in self.records if r.name == record_name]
        if len(lines) > 0:
            return delim.join([l[line_start:].rstrip() for l in lines]).replace("  ", " ")




class TitleSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process HEADER
        header = [r for r in self.records if r.name == "HEADER"]
        if len(header) == 0:
            self.classification, self.date, self.id_code = None, None, None
        else:
            header = header[0]
            self.classification = header[10:50].strip() if header[10:50] != " " * 40 else None
            self.date = header[50:59].strip() if header[50:59] != " " * 9 else None
            self.id_code = header[62:66] if header[62:66] != "    " else None


        #Process OBSLTEs
        obsoletes = [r for r in records if r.name == "OBSLTE"]
        self.obsolete =  (len(obsoletes) > 0)


        #Process TITLEs
        self.title = self._get_record_contents("TITLE")


        #Process SPLITs
        splits = self._get_record_contents("SPLIT")
        self.splits = splits.split() if splits != None else []


        #Process CAVEATs
        caveats = self._get_record_contents("CAVEAT")
        self.caveats = ["%s: %s" % (c[11:15], c[19:79]) for c in caveats] if caveats != None else []


        #Process COMPNDs
        compounds = self._get_record_contents("COMPND", delim="")
        if compounds is None:
            self.compounds = []
        else:
            mol_starts = get_substring_positions(compounds, "MOL_ID")
            mol_starts = sorted(mol_starts)
            molecules = break_into_sections(compounds, mol_starts)
            self.compounds = []
            for molecule in molecules:
                mol = {}
                for identifier in [i for i in molecule.strip().split(";") if i != ""]:
                    key = identifier.split(":")[0].strip()
                    if key not in Pdb.VALID_COMPNDS:
                        pass#raise NotValidPdbError("%s is not a valid idnetifier for COMPND" % key)
                    value = identifier.split(":")[1].strip() if ":" in identifier else ""
                    if "," in value and (key == "CHAIN" or key == "SYNONYM"):
                        value = value.replace(", ",",").split(",")
                    elif value == "YES":
                        value = True
                    elif value == "NO":
                        value = False
                    mol[key] = value
                self.compounds.append(mol)


        #Process SOURCEs
        sources = self._get_record_contents("SOURCE")
        if sources is None:
            self.sources = []
        else:
            mol_starts = get_substring_positions(sources, "MOL_ID")
            mol_starts = sorted(mol_starts)
            molecules = break_into_sections(sources, mol_starts)
            self.sources = []
            for molecule in molecules:
                mol = {}
                for identifier in [i for i in molecule.strip().split(";") if i != ""]:
                    key = identifier.split(":")[0].strip()
                    if key not in Pdb.VALID_SOURCES:
                        pass#raise NotValidPdbError("%s is not a valid idnetifier for SOURCE" % key)
                    value = identifier.split(":")[1].strip() if ":" in identifier else ""
                    if value == "YES":
                        value = True
                    elif value == "NO":
                        value = False
                    mol[key] = value
                self.sources.append(mol)


        #Process KEYWDS
        keywords = self._get_record_contents("KEYWDS")
        if keywords is None:
            self.keywords = None
        else:
            self.keywords = [word.strip() for word in keywords.replace(", ", ",").split(",")]


        #Process EXPDTA
        expdata = self._get_record_contents("EXPDTA", "\n")
        if expdata is None:
            self.expdata = []
        else:
            self.expdata = [d.strip() for d in expdata.split("\n")]
            for datum in self.expdata:
                for keyword in datum.split(";"):
                    if keyword.strip() not in Pdb.VALID_EXPDATA:
                        raise NotValidPdbError("%s is not valid EXPDATA" % keyword)


        #Process NUMMDLs
        model_num = [r for r in self.records if r.name == "NUMMDL"]
        if len(model_num) == 0:
            self.model_num = 1
        else:
            model_num = model_num[0]
            self.model_num = int(model_num[10:14].strip()) if model_num[10:14] != "    " else 1


        #Process MDLTYPs
        #TODO


        #Process AUTHORs
        authors = self._get_record_contents("AUTHOR", "")
        if authors is None:
            self.authors = []
        else:
            self.authors = [a.strip() for a in authors.split(",")]


        #Process REVDATs
        #TODO


        #Process SPRSDEs
        #TODO


        #Process JRNLs
        journals = [r for r in self.records if r.name == "JRNL"]
        if journals is None:
            self.primary_citation = None
        else:
            self.primary_citation = {}

            auths = [line.text for line in journals if line[12:16] == "AUTH"]
            titls = [line.text for line in journals if line[12:16] == "TITL"]
            edits = [line.text for line in journals if line[12:16] == "EDIT"]
            refs = [line.text for line in journals if line[12:15] == "REF"]
            publs = [line.text for line in journals if line[12:16] == "PUBL"]
            refns = [line.text for line in journals if line[12:16] == "REFN"]
            pmids = [line.text for line in journals if line[12:16] == "PMID"]
            dois = [line.text for line in journals if line[12:15] == "DOI"]

            if auths != []:
                self.primary_citation["authors"] = "".join([
                l.strip()[19:] for l in auths
                ]).split(",")
            else:
                self.primary_citation["authors"] = []

            if titls != []:
                self.primary_citation["title"] = " ".join([
                l.strip()[19:].strip() for l in titls
                ])
            else:
                self.primary_citation["authors"] = None

            if edits != []:
                self.primary_citation["editors"] = "".join([
                l.strip()[19:].strip() for l in edits
                ]).split(",")
            else:
                self.primary_citation["editors"] = []

            if refs != []:
                pass
                #TODO

            if publs != []:
                pass
                #TODO

            if refns != []:
                pass
                #TODO

            if pmids != []:
                self.primary_citation["PMID"] = pmids[0][19:].strip()
            else:
                self.primary_citation["PMID"] = None

            if dois != []:
                self.primary_citation["DOI"] = "".join([
                l.strip()[19:] for l in dois
                ])
            else:
                self.primary_citation["DOI"] = None




class PrimaryStructureSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)


        #Process DBREFs (overviews of chains and pointers to reference sequence)
        dbrefs = [r.text for r in self.records if r.name == "DBREF"]
        dbref1 = [r.text for r in self.records if r.name == "DBREF1"]
        dbref2 = [r.text for r in self.records if r.name == "DBREF2"]

        self.dbrefs = []
        for dbref in dbrefs:
            chain = {}
            chain["id"] = dbref[7:11].strip() if dbref[7:11] != "    " else None
            chain["chain"] = dbref[12] if dbref[12] != " " else None
            chain["sequence_start"] = int(dbref[14:18].strip()) if dbref[14:18] != "    " else None
            chain["insertion_start"] = dbref[18] if dbref[18] != " " else None
            chain["sequence_end"] = int(dbref[20:24].strip()) if dbref[20:24] != "    " else None
            chain["insertion_end"] = dbref[24] if dbref[24] != " " else None
            chain["db_name"] = dbref[26:32].strip() if dbref[26:32] != " " * 6 else None
            chain["accession"] = dbref[33:41].strip() if dbref[33:41] != " " * 8 else None
            chain["db_id"] = dbref[42:54].strip() if dbref[42:54] != " " * 12 else None
            chain["db_start"] = int(dbref[55:60].strip()) if dbref[55:60] != "     " else None
            chain["db_insert_start"] = dbref[60] if dbref[60] != " " else None
            chain["db_end"] = int(dbref[62:67]) if dbref[62:67] != "     " else None
            chain["db_insert_end"] = dbref[67] if dbref[67] != " " else None
            self.dbrefs.append(chain)

        assert len(dbref1) == len(dbref2)
        for dbref in zip(dbref1, dbref2):
            chain = {}
            chain["id"] = dbref[0][7:11].strip() if dbref[0][7:11] != "    " else None
            chain["chain"] = dbref[0][12] if dbref[0][12] != " " else None
            chain["sequence_start"] = int(dbref[0][14:18].strip()) if dbref[0][14:18] != "    " else None
            chain["insertion_start"] = dbref[0][18] if dbref[0][18] != " " else None
            chain["sequence_end"] = int(dbref[0][20:24].strip()) if dbref[0][20:24] != "    " else None
            chain["insertion_end"] = dbref[0][24] if dbref[0][24] != " " else None
            chain["db_name"] = dbref[0][26:32].strip() if dbref[0][26:32] != " " * 6 else None
            chain["accession"] = dbref[1][18:40].strip() if dbref[1][18:40] != " " * 22 else None
            chain["db_id"] = dbref[0][48:67].strip() if dbref[0][48:67] != " " * 19 else None
            chain["db_start"] = int(dbref[1][45:55].strip()) if dbref[1][45:55] != " " * 10 else None
            chain["db_insert_start"] = None
            chain["db_end"] = int(dbref[1][57:67]) if dbref[1][57:67] != " " * 10 else None
            chain["db_insert_end"] = None
            self.dbrefs.append(chain)


        #Process SEQADVs (any differences between ref sequence and this sequence)
        seqadvs = [r for r in self.records if r.name == "SEQADV"]
        self.sequence_differences = [{
         "residue": s[12:15].strip() if s[12:15] != "   " else None,
         "chain": s[16] if s[16] != " " else None,
         "sequence_number": int(s[18:22].strip()) if s[18:22] != "    " else None,
         "db_name": s[24:28].strip() if s[24:28] != "    " else None,
         "db_accession": s[29:38].strip() if s[29:38] != " " * 9 else None,
         "db_residue": s[39:42] if s[39:42] != "   " else None,
         "db_sequence_number": int(s[43:48].strip()) if s[43:48] != "     " else None,
         "comment": s[49:70].strip() if s[49:70] != " " * 21 else None
        } for s in seqadvs]


        #Process SEQRESs (the actual sequence)
        seqres = [r.text for r in self.records if r.name == "SEQRES"]
        if seqres is None:
            self.chain_sequences = []
        else:
            chains = list(set([l[11] for l in seqres]))
            self.chain_sequences = []
            for chain in chains:
                lines = [l for l in seqres if l[11] == chain]
                lines = [l[19:70].strip() for l in lines]
                sequence = " ".join(lines)
                self.chain_sequences.append({
                 chain: [r.strip() for r in sequence.split()]
                })


        #Process MODRESs (modified residues)
        modres = [r for r in self.records if r.name == "MODRES"]
        if len(modres) == 0:
            self.modified_residues = []
        else:
            self.modified_residues = []
            for m in modres:
                self.modified_residues.append({
                 "new_residue": m[12:15].strip() if m[12:15] != "   " else None,
                 "chain": m[16] if m[16] != " " else None,
                 "sequence_number": int(m[18:22].strip()) if m[18:22] != "    " else None,
                 "insert_code": m[22] if m[22] != " " else None,
                 "standard_residue": m[24:27].strip() if m[24:27] != "   " else None,
                 "comment": m[29:70].strip() if m[29:70] != " " * 41 else None
                })




class HeterogenSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)


        #Process HETNAMs
        hetnams = [r for r in self.records if r.name == "HETNAM"]
        if len(hetnams) == 0:
            self.hets = []
        else:
            hets = list(set([l[11:14] for l in hetnams]))
            het_dicts = []
            for het in hets:
                name = " ".join(
                 [l[4:].strip() for l in
                  self._get_record_contents("HETNAM", "\n", 11).split("\n")
                   if l[:3] == het]).strip().replace("  ", " ")

                het_dicts.append({"het_id":het, "name":name, "occurences":[], "synonyms":[], "formula":None})
            self.hets = het_dicts

            #Process HETs (Yes it's supposed to be indented - only happens if hetnames are defined)
            hets = [r for r in self.records if r.name == "HET"]
            for het in hets:
                for het_dict in self.hets:
                    if het[7:10].strip() == het_dict["het_id"]:
                        het_dict["occurences"].append({
                         "chain": het[12] if het[12] != " " else None,
                         "sequence_number": int(het[13:17].strip()) if het[13:17] != "   " else None,
                         "insertion_number": het[17] if het[17] != " " else None,
                         "atoms": int(het[20:25].strip()) if het[20:25] != "     " else None,
                         "description": het[30:70].strip() if het[30:70] != " " * 40 else None
                        })

            #Process HETSYNs
            hetsyns = [r.text for r in self.records if r.name == "HETSYN"]
            for het in hetsyns:
                for het_dict in self.hets:
                    if het[11:15].strip() == het_dict["het_id"]:
                        synonyms = [s.strip() for s in het.split(";")]
                        het_dict["synonyms"] += synonyms


        #Process FORMULs
        self.other_molecules = []
        formuls = [r.text for r in self.records if r.name == "FORMUL"]
        mols = list(set([r[12:15].strip() for r in formuls]))

        for mol in mols:
            formula = " ".join(
             [r[19:70].strip() for r in formuls if r[12:15].strip() == mol]
            ).replace("  ", " ")

            if mol in [h["het_id"] for h in self.hets]:
                #This is a het
                for het in self.hets:
                    if het["het_id"] == mol:
                        het["formula"] = formula
            else:
                #This is some random molecule, probably a solvent
                water = [r for r in formuls if r[12:15].strip() == mol][0][18] == "*"
                self.other_molecules.append({
                 "name": mol,
                 "water": water,
                 "formula": formula
                })




class SecondaryStructureSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)


        #Process HELIXs
        helices = [r for r in self.records if r.name == "HELIX"]
        self.helices = [{
         "serial_num": int(l[7:10].strip()) if l[7:10] != "   " else None,
         "helix_id": l[11:14].strip() if l[11:14] != "   " else None,
         "start_residue_name": l[15:18].strip() if l[15:18] != "   " else None,
         "start_residue_chain": l[19] if l[19] != " " else None,
         "start_residue_number": int(l[21:25].strip()) if l[21:25] != "    " else None,
         "start_residue_insert": l[25] if l[25] != " " else None,
         "end_residue_name": l[27:30].strip() if l[27:28] != "   " else None,
         "end_residue_chain": l[31] if l[31] != " " else None,
         "end_residue_number": int(l[33:37].strip()) if l[33:37] != "    " else None,
         "end_residue_insert": l[37] if l[37] != " " else None,
         "helix_class": int(l[38:40].strip()) if l[38:40] != "  " else None,
         "comment": l[40:70].strip() if l[40:70] != " " * 30 else None,
         "length": int(l[71:76].strip()) if l[71:76] != "    " else None
        } for l in helices]


        #Process SHEETs
        sheet_lines = [r for r in self.records if r.name == "SHEET"]
        sheets = list(set([l[11:14].strip() for l in sheet_lines]))
        self.sheets = []

        for sheet in sheets:
            lines = [l for l in sheet_lines if l[11:14].strip() == sheet]
            strands = [{
             "strand_id": l[7:10].strip() if l[7:10] != "   " else None,
             "start_residue_name": l[17:20].strip() if l[17:20] != "   " else None,
             "start_residue_chain": l[21] if l[21] != " " else None,
             "start_residue_number": int(l[22:26].strip()) if l[22:26] != "    " else None,
             "start_residue_insert": l[26] if l[26] != " " else None,
             "end_residue_name": l[28:31].strip() if l[28:31] != "   " else None,
             "end_residue_chain": l[32] if l[32] != " " else None,
             "end_residue_number": int(l[33:37].strip()) if l[33:37] != "    " else None,
             "end_residue_insert": l[37] if l[37] != " " else None,
             "sense": int(l[38:40].strip()) if l[38:40] != "  " else None,
             "reg_cur_atom": l[41:45].strip() if l[41:45] != "    " else None,
             "reg_cur_residue": l[45:48] if l[45:48] != "   " else None,
             "reg_cur_chain": l[49] if l[49] != " " else None,
             "reg_cur_number": int(l[50:54].strip()) if l[50:54] != "    " else None,
             "reg_cur_insert": l[54] if l[54] != " " else None,
             "reg_prev_atom": l[56:60].strip() if l[56:60] != "    " else None,
             "reg_prev_residue": l[60:63] if l[60:63] != "   " else None,
             "reg_prev_chain": l[64] if l[64] != " " else None,
             "reg_prev_number": int(l[65:69].strip()) if l[65:69] != "    " else None,
             "reg_prev_insert": l[69] if l[69] else None
            } for l in lines]
            self.sheets.append({"sheet_id":sheet, "strands":strands})




class ConnectAnnotationSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process SSBONDs
        ssbonds = [r for r in self.records if r.name == "SSBOND"]
        self.ssbonds = [{
         "serial_num": int(s[7:10].strip()) if s[7:10] != "   " else None,
         "residue_1_name": s[11:14].strip() if s[11:14] != "   " else None,
         "residue_1_chain": s[15] if s[15] != " " else None,
         "residue_1_number": int(s[17:21].strip()) if s[17:21] != "    " else None,
         "residue_1_insert": s[21] if s[21] != " " else None,
         "residue_1_symmetry": s[59:65].strip() if s[59:65] != " " * 6 else None,
         "residue_2_name": s[25:28].strip() if s[25:28] != "   " else None,
         "residue_2_chain": s[29] if s[29] != " " else None,
         "residue_2_number": int(s[31:35].strip()) if s[31:35] != "    " else None,
         "residue_2_symmetry": s[66:72].strip() if s[59:65] != " " * 6 else None,
         "residue_2_insert": s[35] if s[35] != " " else None,
         "disulfide_distance": float(s[73:78].strip()) if s[73:78] != "     " else None
        } for s in ssbonds]


        #Process LINKs
        links = [r for r in self.records if r.name == "LINK"]
        self.links = [{
         "residue_1_atom": s[12:16].strip() if s[12:16] != "    " else None,
         "residue_1_name": s[17:20].strip() if s[17:20] != "   " else None,
         "residue_1_chain": s[21] if s[21] != " " else None,
         "residue_1_number": int(s[22:26].strip()) if s[22:26] != "    " else None,
         "residue_1_symmetry": s[59:65].strip() if s[59:65] != " " * 6 else None,
         "residue_1_insert": s[21] if s[21] != " " else None,
         "residue_2_atom": s[42:46].strip() if s[42:46] != "    " else None,
         "residue_2_name": s[47:50].strip() if s[47:50] != "   " else None,
         "residue_2_chain": s[51] if s[51] != " " else None,
         "residue_2_number": int(s[52:56].strip()) if s[52:56] != "    " else None,
         "residue_2_symmetry": s[66:72].strip() if s[59:65] != " " * 6 else None,
         "link_distance": float(s[73:78].strip()) if s[73:78] != "     " else None
        } for s in links]


        #Process CISPEPs
        cispeps = [r for r in self.records if r.name == "CISPEP"]
        self.cispeps = [{
         "serial_num": int(s[7:10].strip()) if s[7:10] != "   " else None,
         "residue_1_name": s[11:14].strip() if s[11:14] != "   " else None,
         "residue_1_chain": s[15] if s[15] != " " else None,
         "residue_1_number": int(s[17:21].strip()) if s[17:21] != "    " else None,
         "residue_1_insert": s[21] if s[21] != " " else None,
         "residue_2_name": s[25:28].strip() if s[25:28] != "   " else None,
         "residue_2_chain": s[29] if s[29] != " " else None,
         "residue_2_number": int(s[31:35].strip()) if s[31:35] != "    " else None,
         "residue_2_insert": s[35] if s[35] != " " else None,
         "modnum": int(s[43:46].strip()) if s[43:46] != "   " else None,
         "angle_measure": float(s[53:59].strip()) if s[53:59] != "     " else None
        } for s in cispeps]




class MiscellaneousSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process SITEs
        sites = [r for r in self.records if r.name == "SITE"]
        site_names = list(set([s[11:14].strip() for s in sites]))
        self.sites = []
        for name in site_names:
            lines = [r for r in sites if r[11:14].strip() == name]
            site = {"name": name, "residues": []}
            for line in lines:
                if line[18:27] != " " * 9:
                    site["residues"].append({
                     "chain": line[22] if line[22] != " " else None,
                     "residue_name": line[18:21].strip() if line[18:21] != "   " else None,
                     "residue_number": int(line[23:27].strip()) if line[23:27] != "    " else None
                    })
                if line[29:38] != " " * 9:
                    site["residues"].append({
                     "chain": line[33] if line[33] != " " else None,
                     "residue_name": line[29:32].strip() if line[29:32] != "   " else None,
                     "residue_number": int(line[34:38].strip()) if line[34:38] != "    " else None
                    })
                if line[40:49] != " " * 9:
                    site["residues"].append({
                     "chain": line[44] if line[44] != " " else None,
                     "residue_name": line[40:43].strip() if line[40:43] != "   " else None,
                     "residue_number": int(line[45:49].strip()) if line[45:49] != "    " else None
                    })
                if line[51:60] != " " * 9:
                    site["residues"].append({
                     "chain": line[55] if line[55] != " " else None,
                     "residue_name": line[51:54].strip() if line[51:54] != "   " else None,
                     "residue_number": int(line[56:60].strip()) if line[56:60] != "    " else None
                    })
            self.sites.append(site)




class CrystalSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process CRYST1
        cryst1 = [r for r in self.records if r.name == "CRYST1"]
        if len(cryst1) == 0:
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.s_group, self.z = (
             None, None, None, None, None, None, None, None)
        else:
            cryst1 = cryst1[0]
            self.a = float(cryst1[6:15].strip()) if cryst1[6:15] != " " * 9 else None
            self.b = float(cryst1[15:24].strip()) if cryst1[15:24] != " " * 9 else None
            self.c = float(cryst1[24:33].strip()) if cryst1[24:33] != " " * 9 else None
            self.alpha = float(cryst1[33:40].strip()) if cryst1[33:40] != " " * 7 else None
            self.beta = float(cryst1[40:47].strip()) if cryst1[40:47] != " " * 7 else None
            self.gamma = float(cryst1[47:54].strip()) if cryst1[47:54] != " " * 7 else None
            self.s_group = cryst1[55:66].strip() if cryst1[55:66] != " " * 11 else None
            self.z = int(cryst1[66:70].strip()) if cryst1[66:70] != " " * 4 else None


        #Process ORIGXns
        origx1 = [r for r in records if r.name == "ORIGX1"]
        origx2 = [r for r in records if r.name == "ORIGX2"]
        origx3 = [r for r in records if r.name == "ORIGX3"]

        if len(origx1) == 0:
            self.o11, self.o12, self.o13, self.t1 = None, None, None, None
        else:
            origx1 = origx1[0]
            self.o11 = float(origx1[10:20].strip()) if origx1[10:20] != " " * 10 else None
            self.o12 = float(origx1[20:30].strip()) if origx1[20:30] != " " * 10 else None
            self.o13 = float(origx1[30:40].strip()) if origx1[30:40] != " " * 10 else None
            self.t1 = float(origx1[45:55].strip()) if origx1[45:55] != " " * 10 else None

        if len(origx2) == 0:
            self.o21, self.o22, self.o23, self.t2 = None, None, None, None
        else:
            origx2 = origx2[0]
            self.o21 = float(origx2[10:20].strip()) if origx2[10:20] != " " * 10 else None
            self.o22 = float(origx2[20:30].strip()) if origx2[20:30] != " " * 10 else None
            self.o23 = float(origx2[30:40].strip()) if origx2[30:40] != " " * 10 else None
            self.t2 = float(origx2[45:55].strip()) if origx2[45:55] != " " * 10 else None

        if len(origx3) == 0:
            self.o31, self.o32, self.o33, self.t3 = None, None, None, None
        else:
            origx3 = origx3[0]
            self.o31 = float(origx3[10:20].strip()) if origx3[10:20] != " " * 10 else None
            self.o32 = float(origx3[20:30].strip()) if origx3[20:30] != " " * 10 else None
            self.o33 = float(origx3[30:40].strip()) if origx3[30:40] != " " * 10 else None
            self.t3 = float(origx3[45:55].strip()) if origx3[45:55] != " " * 10 else None


        #Process SCALEns
        scale1 = [r for r in records if r.name == "SCALE1"]
        scale2 = [r for r in records if r.name == "SCALE2"]
        scale3 = [r for r in records if r.name == "SCALE3"]

        if len(scale1) == 0:
            self.s11, self.s12, self.s13, self.u1 = None, None, None, None
        else:
            scale1 = scale1[0]
            self.s11 = float(scale1[10:20].strip()) if scale1[10:20] != " " * 10 else None
            self.s12 = float(scale1[20:30].strip()) if scale1[20:30] != " " * 10 else None
            self.s13 = float(scale1[30:40].strip()) if scale1[30:40] != " " * 10 else None
            self.u1 = float(scale1[45:55].strip()) if scale1[45:55] != " " * 10 else None

        if len(scale2) == 0:
            self.s21, self.s22, self.s23, self.u2 = None, None, None, None
        else:
            scale2 = scale2[0]
            self.s21 = float(scale2[10:20].strip()) if scale2[10:20] != " " * 10 else None
            self.s22 = float(scale2[20:30].strip()) if scale2[20:30] != " " * 10 else None
            self.s23 = float(scale2[30:40].strip()) if scale2[30:40] != " " * 10 else None
            self.u2 = float(scale2[45:55].strip()) if scale2[45:55] != " " * 10 else None

        if len(scale3) == 0:
            self.s31, self.s32, self.s33, self.u3 = None, None, None, None
        else:
            scale3 = scale3[0]
            self.s31 = float(scale3[10:20].strip()) if scale3[10:20] != " " * 10 else None
            self.s32 = float(scale3[20:30].strip()) if scale3[20:30] != " " * 10 else None
            self.s33 = float(scale3[30:40].strip()) if scale3[30:40] != " " * 10 else None
            self.u3 = float(scale3[45:55].strip()) if scale3[45:55] != " " * 10 else None


        #Process MTRIXns
        mtrix1 = [r for r in records if r.name == "MTRIX1"]
        mtrix2 = [r for r in records if r.name == "MTRIX2"]
        mtrix3 = [r for r in records if r.name == "MTRIX3"]

        if len(mtrix1) == 0:
            self.m11, self.m12, self.m13, self.v1 = None, None, None, None
        else:
            mtrix1 = mtrix1[0]
            self.m11 = float(mtrix1[10:20].strip()) if mtrix1[10:20] != " " * 10 else None
            self.m12 = float(mtrix1[20:30].strip()) if mtrix1[20:30] != " " * 10 else None
            self.m13 = float(mtrix1[30:40].strip()) if mtrix1[30:40] != " " * 10 else None
            self.v1 = float(mtrix1[45:55].strip()) if mtrix1[45:55] != " " * 10 else None

        if len(mtrix2) == 0:
            self.m21, self.m22, self.m23, self.v2 = None, None, None, None
        else:
            mtrix2 = mtrix2[0]
            self.m21 = float(mtrix2[10:20].strip()) if mtrix2[10:20] != " " * 10 else None
            self.m22 = float(mtrix2[20:30].strip()) if mtrix2[20:30] != " " * 10 else None
            self.m23 = float(mtrix2[30:40].strip()) if mtrix2[30:40] != " " * 10 else None
            self.v2 = float(mtrix2[45:55].strip()) if mtrix2[45:55] != " " * 10 else None

        if len(mtrix3) == 0:
            self.m31, self.m32, self.m33, self.v3 = None, None, None, None
        else:
            mtrix3 = mtrix3[0]
            self.m31 = float(mtrix3[10:20].strip()) if mtrix3[10:20] != " " * 10 else None
            self.m32 = float(mtrix3[20:30].strip()) if mtrix3[20:30] != " " * 10 else None
            self.m33 = float(mtrix3[30:40].strip()) if mtrix3[30:40] != " " * 10 else None
            self.v3 = float(mtrix3[45:55].strip()) if mtrix3[45:55] != " " * 10 else None




class CoordinateSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process MODELs
        models = [r for r in self.records if r.name == "MODEL"]

        if len(models) != 0:
            #There are multiple models here
            models = [x[1:-1] for x in
             break_into_sections(self.records, [self.records.index(m) for m in models])]
        else:
            #There is just one model here
            models = [self.records]
        #Each object in models should now be an array of records

        self.models = []
        for model_lines in models:
            model = {}

            #Process ATOMs
            atoms = [r for r in model_lines if r.name == "ATOM"]
            model["atoms"] = [{
             "serial_num": int(a[6:11].strip()) if a[6:11] != "     " else None,
             "atom_name": a[12:16].strip() if a[12:16] != "    " else None,
             "alt_loc": a[16] if a[16] != " " else None,
             "residue_name": a[17:20].strip() if a[17:20] != "   " else None,
             "chain": a[21] if a[21] != " " else None,
             "residue_number": int(a[22:26].strip()) if a[22:26] != "    " else None,
             "residue_insert": a[26] if a[26] != " " else None,
             "x": float(a[30:38].strip()) if a[30:38] != " " * 8 else None,
             "y": float(a[38:46].strip()) if a[38:46] != " " * 8 else None,
             "z": float(a[46:54].strip()) if a[46:54] != " " * 8 else None,
             "occupancy": float(a[54:60].strip()) if a[54:60] != " " * 6 else None,
             "temp_factor": float(a[60:66].strip()) if a[60:66] != " " * 6 else None,
             "element": a[76:78].strip() if a[76:78] != "  " else None,
             "charge": a[78:80].strip() if a[78:80] != "  " else None
            } for a in atoms]


            #Process ANISOU
            anisou = [r for r in model_lines if r.name == "ANISOU"]
            model["anisou"] = [{
             "serial_num": int(a[6:11].strip()) if a[6:11] != "     " else None,
             "atom_name": a[12:16].strip() if a[12:16] != "    " else None,
             "alt_loc": a[16] if a[16] != " " else None,
             "residue_name": a[17:20].strip() if a[17:20] != "   " else None,
             "chain": a[21] if a[21] != " " else None,
             "residue_number": int(a[22:26].strip()) if a[22:26] != "    " else None,
             "residue_insert": a[26] if a[26] != " " else None,
             "u11": int(a[28:35].strip()) if a[28:35] != " " * 7 else None,
             "u22": int(a[35:42].strip()) if a[35:42] != " " * 7 else None,
             "u33": int(a[42:49].strip()) if a[42:49] != " " * 7 else None,
             "u12": int(a[49:56].strip()) if a[49:56] != " " * 7 else None,
             "u13": int(a[56:63].strip()) if a[56:63] != " " * 7 else None,
             "u23": int(a[63:70].strip()) if a[63:70] != " " * 7 else None,
             "element": a[76:78].strip() if a[76:78] != "  " else None,
             "charge": a[78:80].strip() if a[78:80] != "  " else None
            } for a in anisou]


            #Process HETATMs
            hetatoms = [r for r in model_lines if r.name == "HETATM"]
            model["hetero_atoms"] = [{
             "serial_num": int(a[6:11].strip()) if a[6:11] != "     " else None,
             "atom_name": a[12:16].strip(),
             "alt_loc": a[16] if a[16] != " " else None,
             "residue_name": a[17:20] if a[17:20] != "   " else None,
             "chain": a[21] if a[21] != " " else None,
             "residue_number": int(a[22:26].strip()) if a[22:26] != "    " else None,
             "x": float(a[30:38].strip()),
             "y": float(a[38:46].strip()),
             "z": float(a[46:54].strip()),
             "occupancy": float(a[54:60].strip()) if a[54:60] != " " * 6 else None,
             "temp_factor": float(a[60:66].strip()) if a[60:66] != " " * 6 else None,
             "element": a[76:78].strip() if a[76:78] != "  " else None,
             "charge": a[78:80].strip() if a[78:80] != "  " else None
            } for a in hetatoms]


            self.models.append(model)




class ConnectSection(Section):

    def __init__(self, records):
        Section.__init__(self, records)

        #Process CONECTs
        atoms = list(set([int(r[6:11].strip()) for r in self.records]))

        self.atoms = []
        for atom_id in atoms:
            atom = {"atom_id": atom_id, "bonded_atoms": []}
            strings = [r[11:31] for r in self.records if r[6:11].strip() == str(atom_id)]
            for string in strings:
                if string[:5] != "     ":
                    atom["bonded_atoms"].append(int(string[:5].strip()))
                if string[5:10] != "     ":
                    atom["bonded_atoms"].append(int(string[5:10].strip()))
                if string[10:15] != "     ":
                    atom["bonded_atoms"].append(int(string[10:15].strip()))
                if string[15:20] != "     ":
                    atom["bonded_atoms"].append(int(string[15:20].strip()))
            self.atoms.append(atom)




class Model:
    """This class represents the structure contained within a PDB file"""

    def __init__(self, model_dict, sites):
        #Get chains
        chains = list(set([a["chain"] for a in model_dict["atoms"]]))
        self.chains = []

        for chain in chains:
            atoms = [a for a in model_dict["atoms"] if a["chain"] == chain]
            anisous = [a for a in model_dict["anisou"] if a["chain"] == chain]
            self.chains.append(Chain(atoms, anisous))


        #Get heteroatoms
        hets = list(set([(h["residue_number"], h["chain"]) for h in model_dict["hetero_atoms"]]))
        self.hets = []

        for het in hets:
            atoms = [a for a in model_dict["hetero_atoms"]
             if a["residue_number"] == het[0] and a["chain"] == het[1]]
            anisous = [a for a in model_dict["anisou"]
             if a["residue_number"] == het[0] and a["chain"] == het[1]]
            self.hets.append(Het(atoms, anisous))


        #Process chains
        self.mass = 0
        self.residues = []
        self.atoms = []

        for chain in self.chains:
            chain.model = self
            self.residues += chain.residues
            self.atoms += chain.atoms
            self.mass += chain.mass


        #Process hets
        for het in self.hets:
            het.model = self
            self.mass += het.mass
            self.atoms += het.atoms


        #Process residues
        for residue in self.residues:
            residue.model = self


        #Process atoms
        for atom in self.atoms:
            atom.model = self


        #Add in sites
        self.sites = []
        for site in sites:
            residues = []
            for residue in site["residues"]:
                for chain in self.chains:
                    if chain.name == residue["chain"]:
                        for obj_residue in chain.residues:
                            if obj_residue.number == residue["residue_number"]:
                                residues.append(obj_residue)
            if len(residues) > 0:
                self.sites.append(Site(site["name"], residues))





class Chain:
    """This class represents a PDB chain"""

    def __init__(self, atoms, anisous):
        self.name = atoms[0]["chain"]

        #Get residues
        residues = list(set([a["residue_number"] for a in atoms]))
        self.residues = []

        for residue in residues:
            res_atoms = [a for a in atoms if a["residue_number"] == residue]
            res_anisous = [a for a in anisous if a["residue_number"] == residue]
            self.residues.append(Residue(res_atoms, res_anisous))

        #Process residues
        self.mass = 0
        self.atoms = []
        for residue in self.residues:
            residue.chain = self
            self.atoms += residue.atoms
            self.mass += residue.mass


        #Process atoms
        for atom in self.atoms:
            atom.chain = self




class Residue:
    """This class represents a PDB residue"""

    def __init__(self, atoms, anisous):
        self.name = atoms[0]["residue_name"]
        self.number = atoms[0]["residue_number"]

        #Get atoms
        unique_atoms =  list(set([a["serial_num"] for a in atoms]))
        self.atoms = []

        for atom in unique_atoms:
            self.atoms.append(Atom(
             [a for a in atoms if a["serial_num"] == atom][0],
             [a for a in anisous if a["serial_num"] == atom]
            ))

        #Process atoms
        self.mass = 0
        for atom in self.atoms:
            atom.residue = self
            self.mass += PERIODIC_TABLE[atom.element.upper()]




class Het(Residue):

    def __init__(self, atoms, anisous):
        Residue.__init__(self, atoms, anisous)




class Atom:
    """This class represents an atom"""

    def __init__(self, atom, anisous):
        self.number = atom["serial_num"]
        self.name = atom["atom_name"]

        self.x = atom["x"]
        self.y = atom["y"]
        self.z = atom["z"]

        self.occupancy = atom["occupancy"]
        self.temp_factor = atom["temp_factor"]

        self.element = atom["element"]
        self.charge = atom["charge"]

        self.residue = None


    def distance_to(self, other_atom):
        x_sum = math.pow((other_atom.x - self.x), 2)
        y_sum = math.pow((other_atom.y - self.y), 2)
        z_sum = math.pow((other_atom.z - self.z), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)

        return distance




class Site:

    def __init__(self, name, residues):
        self.name = name
        self.residues = residues



    def get_unbroken_site(self):
        #Are all the residues on the same chain?
        chain_of_first = self.residues[0].chain
        for residue in self.residues:
            if residue.chain is not chain_of_first:
                #No, this is impossible
                return None


        #They are all on the same chain, get the residues needed to make unbroken
        residues_needed = range(
         min([r.number for r in self.residues]),
         max([r.number for r in self.residues]) + 1
        )
        residues = []
        for number in residues_needed:
            for residue in chain_of_first.residues:
                if residue.number == number:
                    residues.append(residue)

        return Site(self.name + "_unbroken", residues)


    def get_pymol_selector_string(self):
        s = []

        for residue in [r for r in self.residues if type(r) != type(Het)]:
            s.append(
             "(resi %i & chain %s)" % (residue.number, residue.chain.name)
            )

        return " | ".join(s)



    def __repr__(self):
        s = ", ".join([
         "%i (%s)" % (r.number, r.name) for r in self.residues
        ])
        return s



class Pdb:
    VALID_CHARS = """'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ12
    34567890`-=[]\;',./~!@#$%^&*()_+{}|:"<>?'"""

    VALID_RECORDS = ["CRYST1", "END", "HEADER", "NUMMDL", "MASTER", "ORIGX[1-3]", "SCALE[1-3]",
     "AUTHOR", "CAVEAT", "COMPND", "EXPDTA", "MDLTYP", "KEYWDS", "OBSLTE", "SOURCE", "SPLIT", "SPRSDE", "TITLE",
      "ANISOU", "ATOM", "CISPEP", "CONECT", "DBREF", "HELIX", "HET", "HETATM", "LINK", "MODRES", "MTRIX[1-3]", "REVDAT", "SEQADV", "SHEET", "SSBOND",
       "FORMUL", "HETNAM", "HETSYN", "SEQRES", "SITE",
        "ENDMDL", "MODEL", "TER",
         "JRNL", "REMARK",
          "HYDBND", "SLTBRG", "TURN"]

    TITLE_RECORDS = ["HEADER", "OBSLTE", "TITLE", "SPLIT", "CAVEAT", "COMPND",
     "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL"]

    VALID_COMPNDS = ["MOL_ID", "MOLECULE", "CHAIN", "FRAGMENT", "SYNONYM", "EC",
     "ENGINEERED", "MUTATION", "OTHER_DETAILS"]

    VALID_SOURCES = ["MOL_ID", "SYNTHETIC", "FRAGMENT", "ORGANISM_SCIENTIFIC",
     "ORGANISM_COMMON", "ORGANISM_TAXID", "STRAIN", "VARIANT", "CELL_LINE", "ATCC",
      "ORGAN", "TISSUE", "CELL", "ORGANELLE", "SECRETION", "CELLULAR_LOCATION",
       "PLASMID", "GENE", "EXPRESSION_SYSTEM", "EXPRESSION_SYSTEM_COMMON",
        "EXPRESSION_SYSTEM_TAXID", "EXPRESSION_SYSTEM_STRAIN", "EXPRESSION_SYSTEM_VARIANT",
         "EXPRESSION_SYSTEM_CELL_LINE", "EXPRESSION_SYSTEM_ATCC_NUMBER",
          "EXPRESSION_SYSTEM_ORGAN", "EXPRESSION_SYSTEM_TISSUE", "EXPRESSION_SYSTEM_CELL",
           "EXPRESSION_SYSTEM_ORGANELLE", "EXPRESSION_SYSTEM_CELLULAR_LOCATION",
            "EXPRESSION_SYSTEM_VECTOR_TYPE", "EXPRESSION_SYSTEM_VECTOR",
             "EXPRESSION_SYSTEM_PLASMID", "EXPRESSION_SYSTEM_GENE", "OTHER_DETAILS"]

    VALID_EXPDATA = ["X-RAY DIFFRACTION", "FIBER DIFFRACTION", "NEUTRON DIFFRACTION",
     "ELECTRON CRYSTALLOGRAPHY", "ELECTRON MICROSCOPY", "SOLID-STATE NMR",
      "SOLUTION NMR", "SOLUTION SCATTERING"]

    PRIMARY_STRUCTURE_RECORDS = ["DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES"]

    HETEROGEN_RECORDS = ["HET", "HETNAM", "HETSYN", "FORMUL"]

    SECONDARY_STRUCTURE_RECORDS = ["HELIX", "SHEET"]

    CONNECTIVITY_ANNOTATION_RECORDS = ["SSBOND", "LINK", "CISPEP"]

    MISCELLANEOUS_RECORDS = ["SITE"]

    CRYSTAL_RECORDS = ["CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2",
     "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3"]

    COORDINATE_RECORDS = ["MODEL", "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL"]

    CONNECT_RECORDS = ["CONECT"]


    def __init__(self, f, narrate=False):

        #Open file contents
        try:
            contents = f.read()
            contents.replace("\r\n", "\n")
        except UnicodeDecodeError:
            raise NotValidPdbError("This file is not in plain text")


        #Verify that the file characters are valid
        for char in contents:
            if char not in Pdb.VALID_CHARS:
                raise NotValidPdbError("Contains invalid character: " + char)


        #Get lines of file
        self.lines = [line for line in contents.split("\n") if line.strip() != ""]
        if narrate:
            print("There are %i lines in this file." % len(self.lines))
        for line in self.lines[:-1]:
            if len(line) != 80:
                raise NotValidPdbError("Line %i is not 80 characters long" %
                 (self.lines.index(line) + 1))


        #Turn into records
        line_num = 1
        self.records = []
        for line in self.lines:
            self.records.append(Record(line, line_num))
            line_num += 1


        #Get title section
        self.title = TitleSection(
         [r for r in self.records if r.name in Pdb.TITLE_RECORDS])


        #Get primary structure section
        self.primary_structure = PrimaryStructureSection(
         [r for r in self.records if r.name in Pdb.PRIMARY_STRUCTURE_RECORDS])


        #Get heterogen section
        self.heterogen = HeterogenSection(
         [r for r in self.records if r.name in Pdb.HETEROGEN_RECORDS])


        #Get heterogen section
        self.secondary_structure = SecondaryStructureSection(
         [r for r in self.records if r.name in Pdb.SECONDARY_STRUCTURE_RECORDS])


        #Get connectivity annotation section
        self.connect_annotation = ConnectAnnotationSection(
         [r for r in self.records if r.name in Pdb.CONNECTIVITY_ANNOTATION_RECORDS])


        #Get miscellaneous annotation section
        self.miscellaneous = MiscellaneousSection(
         [r for r in self.records if r.name in Pdb.MISCELLANEOUS_RECORDS])


        #Get crystal annotation section
        self.crystal = CrystalSection(
         [r for r in self.records if r.name in Pdb.CRYSTAL_RECORDS])


        #Get coordinate section
        self.coordinate = CoordinateSection(
         [r for r in self.records if r.name in Pdb.COORDINATE_RECORDS])


        #Get connect section
        self.connectivity = ConnectSection(
         [r for r in self.records if r.name in Pdb.CONNECT_RECORDS])


        #Perform verification
        #TODO


        #Get objects
        self.models = [Model(m, self.miscellaneous.sites) for m in self.coordinate.models]
        self.model = self.models[0]


        #Report back on file contents?
        if narrate:
            print("\nMacromolecules contained:")
            if len(self.title.compounds) == 0:
                print("\tNone reported")
            else:
                for compound in self.title.compounds:
                    if "MOLECULE" in compound.keys():
                        print("\t%s" % compound["MOLECULE"])
                    else:
                        print("\t<UNAMED>")
                    if "CHAIN" in compound.keys():
                        print("\t\tChains: %s" % ", ".join(compound["CHAIN"]))
                    else:
                        print("\t\tChains: <None specified>")
                    print("")

            print("\nChains contained:")
            if len(self.primary_structure.dbrefs) == 0:
                print("\tNone reported")
            else:
                for chain in self.primary_structure.dbrefs:
                    print("\t%s (%i residues)" %
                     (chain["chain"], (chain["sequence_end"] - chain["sequence_start"] + 1)))
                    print("\t%i to %i in %s (%s)" %
                     (chain["sequence_start"], chain["sequence_end"], chain["accession"], chain["db_name"]))
                    print("\t%i to %i in db numbering" % (chain["db_start"], chain["db_end"]))
                    print("")

            print("\nNon-Standard Residues (hets) contained:")
            if len(self.heterogen.hets) == 0:
                print("\tNone reported")
            else:
                for het in self.heterogen.hets:
                    print("\t%s (%s)" % (het["het_id"], het["name"]))
                    print("\tFormula: %s" % het["formula"])
                    print("")

            print("\nOther molecules contained:")
            if len(self.heterogen.other_molecules) == 0:
                print("\tNone reported")
            else:
                for mol in self.heterogen.other_molecules:
                    print("\t%s: %s" % (mol["name"], mol["formula"]))
                    print("\tWater: %s" % mol["water"])
                    print("")

            print("\nAnnotated sites:")
            if len(self.miscellaneous.sites) == 0:
                print("\tNone reported")
            else:
                for site in self.miscellaneous.sites:
                    print("\t%s: %i residues" % (site["name"], len(site["residues"])))
                    print("")
