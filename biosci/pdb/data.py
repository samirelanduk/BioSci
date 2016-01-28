from .exceptions import *
import datetime

class PdbDataStructure:
    """A processed PdbFile, mostly in dictionary form"""

    def __init__(self, pdb_file):
        self.file = pdb_file

        #Create the sections
        self.title = TitleSection(self.file)
        self.primary_structure = PrimaryStructureSection(self.file)
        self.heterogen = HeterogenSection(self.file)
        self.secondary_structure = SecondaryStructureSection(self.file)
        self.connectivity_annotation = ConnectivityAnnotationSection(self.file)
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



class TitleSection(PdbSection):

    RECORD_NAMES = ("HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND",
     "SOURCE", "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT",
      "SPRSDE", "JRNL", "REMARK")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process HEADER
        headers = self.get_records_by_name("HEADER")
        if headers:
            header = headers[0]
            self.classification = header[10:50].strip() if header[10:50].strip() else None
            self.date = datetime.datetime.strptime(header[50:59].strip(), "%d-%b-%y").date() if header[50:59].strip() else None
            self.code = header[62:66] if header[62:66].strip() else None
        else:
            self.classification, self.date, self.code = None, None, None

        #Process OBSLTE
        obsolete = self.get_records_by_name("OBSLTE")
        self.obsolete = bool(obsolete)

        #Process TITLE
        title = self.get_records_by_name("TITLE")
        self.title = " ".join([r[10:].strip() for r in title]).strip().replace("  ", " ").strip()

        #Process SPLIT
        splits = self.get_records_by_name("SPLIT")
        self.split_codes = []
        for r in splits:
            for i in range(13):
                code = r[((i+2) * 5) + 1: ((i+2) * 5) + 5]
                if code.strip(): self.split_codes.append(code)

        #Process CAVEATs
        caveat = self.get_records_by_name("CAVEAT")
        self.caveat = " ".join([r[19:].strip() for r in caveat]).strip().replace("  ", " ").strip()

        #Process COMPND
        compnds = self.get_records_by_name("COMPND")
        compound_text = " ".join([r[10:].strip() for r in compnds]).replace("  ", " ").replace("; ", ";").replace(": ", ":")
        pairs = [p.split(":") for p in compound_text.split(";")]
        for offset in range(1, len(pairs))[::-1]:
            if len(pairs[offset]) == 1:
                try:
                    pairs[offset-1][1] += ";" + pairs[offset][0]
                except IndexError:
                    pairs[offset-1][0] += ";" + pairs[offset][0]
        pairs = [p for p in pairs if len(pairs) == 2]
        self.compounds = []
        current_compound = {}
        for pair in pairs:
            if pair[0] == "MOL_ID":
                self.compounds.append(current_compound)
                current_compound = {}
            current_compound[pair[0]] = pair[1]
        self.compounds.append(current_compound)
        self.compounds = [c for c in self.compounds if c]
        for compound in self.compounds:
            for key in compound:
                try:
                    compound[key] = int(compound[key])
                except ValueError:
                    pass
                if compound[key] == "YES":
                    compound[key] = True
                if compound[key] == "NO":
                    compound[key] = False
                if key == "CHAIN":
                    compound[key] = compound[key].replace(" ", "").split(",")

        #Process SOURCE
        sources = self.get_records_by_name("SOURCE")
        source_text = " ".join([r[10:].strip() for r in sources]).replace("  ", " ").replace("; ", ";").replace(": ", ":")
        pairs = [p.split(":") for p in source_text.split(";")]
        for offset in range(1, len(pairs))[::-1]:
            if len(pairs[offset]) == 1:
                try:
                    pairs[offset-1][1] += ";" + pairs[offset][0]
                except IndexError:
                    pairs[offset-1][0] += ";" + pairs[offset][0]
        pairs = [p for p in pairs if len(pairs) == 2]
        self.sources = []
        current_source = {}
        for pair in pairs:
            if pair[0] == "MOL_ID":
                self.sources.append(current_source)
                current_source = {}
            current_source[pair[0]] = pair[1]
        self.sources.append(current_source)
        self.sources = [c for c in self.sources if c]
        for source in self.sources:
            for key in source:
                try:
                    source[key] = int(source[key])
                except ValueError:
                    pass
                if source[key] == "YES":
                    source[key] = True
                if source[key] == "NO":
                    source[key] = False

        #Process KEYWDS
        keywords = self.get_records_by_name("KEYWDS")
        keyword_text = " ".join([r[10:].strip() for r in keywords]).replace("  ", " ").replace(", ", ",")
        self.keywords = keyword_text.split(",")

        #Process EXPDTA
        expdata = self.get_records_by_name("EXPDTA")
        self.experimental_techniques = []
        for r in expdata:
            self.experimental_techniques.append(r[10:].strip().replace("; ", ";").split(";"))

        #Process NUMMDL
        nummdl = self.get_records_by_name("NUMMDL")
        if nummdl:
            self.model_number = int(nummdl[0][10:14].strip()) if nummdl[0][10:14].strip() else 1
        else:
            self.model_number = 1

        #Process MDLTYP
        mdltyp = self.get_records_by_name("MDLTYP")
        self.model_info = []
        current_info = []
        for r in mdltyp:
            if not r[8:10].strip() or int(r[8:10].strip()) == 1:
                self.model_info.append(current_info)
                current_info = []
            current_info.append(r[10:].strip().replace("; ", ";").split(";"))
        self.model_info = [i for i in self.model_info if i]

        #Process AUTHOR
        authors = self.get_records_by_name("AUTHOR")
        author_text = " ".join([r[10:].strip() for r in authors]).replace("  ", " ").replace(", ", ",")
        self.authors = author_text.split(",")

        #Process REVDAT
        revdats = self.get_records_by_name("REVDAT")
        self.modifications = []
        current_modification = {"details": []}
        for r in revdats:
            if not r[10:12].strip() or int(r[10:12].strip()) == 1:
                self.modifications.append(current_modification)
                current_modification = {"details": []}
                current_modification["date"] = datetime.datetime.strptime(r[13:22].strip(), "%d-%b-%y").date()
                current_modification["code"] = r[23:27] if r[23:27].strip() else None
                current_modification["type"] = int(r[31]) if r[31] else 1
            if r[39:45].strip():
                current_modification["details"].append(r[39:45].strip())
            if r[46:52].strip():
                current_modification["details"].append(r[46:52].strip())
            if r[53:59].strip():
                current_modification["details"].append(r[53:59].strip())
            if r[60:66].strip():
                current_modification["details"].append(r[60:66].strip())
        self.modifications = [m for m in self.modifications if "date" in m.keys()]

        #Process SPRSDE
        sprsde = self.get_records_by_name("SPRSDE")
        self.superceded = []
        current_supercede = {"olds": []}
        for r in sprsde:
            if not r[8:10].strip() or int(r[10:12].strip()) == 1:
                self.superceded.append(current_supercede)
                current_supercede = {"olds": []}
                current_supercede["date"] = datetime.datetime.strptime(r[11:20].strip(), "%d-%b-%y").date()
                current_supercede["code"] = r[21:25] if r[23:27].strip() else None
            for s in range(9):
                if r[((s + 6) * 5) + 1: ((s + 6) * 5) + 5].strip():
                    current_supercede["olds"].append(r[((s + 6) * 5) + 1: ((s + 6) * 5) + 5].strip())
        self.superceded = [s for s in self.superceded if "date" in s.keys()]

        #Process JRNL
        jrnl = self.get_records_by_name("JRNL")
        if jrnl:
            self.journal = {}
            author_text = " ".join([r[19:].strip() for r in jrnl if r[12:16] == "AUTH"]
             ).replace("  ", " ").replace(", ", ",")
            self.journal["authors"] = author_text.split(",")
            self.journal["title"] = " ".join([r[19:].strip() for r in jrnl if r[12:16] == "TITL"]
             ).replace("  ", " ")
            editor_text = " ".join([r[19:].strip() for r in jrnl if r[12:16] == "EDIT"]
             ).replace("  ", " ").replace(", ", ",")
            self.journal["editors"] = [e for e in editor_text.split(",") if e]
            ref = [r for r in jrnl if r[12:16].strip() == "REF"]
            if ref:
                self.journal["reference"] = ref[0][19:].strip()
                if "TO BE PUBLISHED" not in self.journal["reference"]:
                    self.journal["reference"] = {}
                    self.journal["reference"]["volume"] = ref[0][51:55].strip() if ref[0][51:55].strip() else None
                    self.journal["reference"]["page"] = ref[0][56:61].strip() if ref[0][56:61].strip() else None
                    self.journal["reference"]["year"] = int(ref[0][62:66].strip()) if ref[0][62:66].strip() else None
                    self.journal["reference"]["publication"] = " ".join([r[19:46].strip() for r in ref]).replace("  ", " ")
            else:
                self.journal["reference"] = None
            publ = [r for r in jrnl if r[12:16] == "PUBL"]
            if publ:
                self.journal["publisher"] = " ".join([r[19:].strip() for r in publ]
                 ).replace("  ", " ")
            else:
                self.journal["publisher"] = None
            refn = [r for r in jrnl if r[12:16] == "REFN"]
            if refn:
                self.journal["reference"] = {}
                self.journal["reference"]["type"] = refn[0][12:16]
                self.journal["reference"]["serial"] = refn[0][40:65].strip()
            else:
                self.journal["reference"] = None
            pmid = [r for r in jrnl if r[12:16] == "PMID"]
            self.journal["pmid"] = pmid[0][19:].strip() if pmid else None
            doi = [r for r in jrnl if r[12:16].strip() == "DOI"]
            self.journal["doi"] = doi[0][19:].strip() if doi else None
        else:
            self.journal = None

        #Process REMARKs
        remarks = self.get_records_by_name("REMARK")
        remark_nums = list(set([int(r[7:10].strip()) for r in remarks]))
        remark_nums.sort()
        self.remarks = []
        for num in remark_nums:
            self.remarks.append({
             "num": num,
             "content": "\n".join([r[11:].rstrip() for r in remarks
              if int(r[7:10].strip()) == num and r[11:].strip()])
            })



class PrimaryStructureSection(PdbSection):

    RECORD_NAMES = ("DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process DBREF
        dbrefs = self.get_records_by_name("DBREF")
        dbrefs_long = zip(self.get_records_by_name("DBREF1"),
         self.get_records_by_name("DBREF2"))
        self.db_refs = [{
         "code": dbref[7:11].strip() if dbref[7:11].strip() else None,
         "chain": dbref[12] if dbref[12].strip() else None,
         "seq_begin": int(dbref[14:18].strip()) if dbref[14:18].strip() else None,
         "insert_begin": dbref[18] if dbref[18].strip() else None,
         "seq_end": int(dbref[20:24].strip()) if dbref[20:24].strip() else None,
         "insert_end": dbref[24] if dbref[24].strip() else None,
         "database": dbref[26:32].strip() if dbref[26:32].strip() else None,
         "accession": dbref[33:41].strip() if dbref[33:41].strip() else None,
         "id": dbref[42:54].strip() if dbref[42:54].strip() else None,
         "db_seq_begin": int(dbref[55:60].strip()) if dbref[55:60].strip() else None,
         "db_insert_begin": dbref[60] if dbref[60].strip() else None,
         "db_seq_end": int(dbref[62:67].strip()) if dbref[62:67].strip() else None,
         "db_insert_end": dbref[67] if dbref[67].strip() else None
        } for dbref in dbrefs]
        for dbref, dbref2 in dbrefs_long:
            self.dbrefs += {
             "code": dbref[7:11].strip() if dbref[7:11] else None,
             "chain": dbref[12] if dbref[12].strip() else None,
             "seq_begin": int(dbref[14:18].strip()) if dbref[14:18].strip() else None,
             "insert_begin": dbref[18] if dbref[18].strip() else None,
             "seq_end": int(dbref[20:24].strip()) if dbref[20:24].strip() else None,
             "insert_end": dbref[24] if dbref[24].strip() else None,
             "database": dbref[26:32].strip() if dbref[26:32].strip() else None,
             "accession": dbref2[18:40].strip() if dbref2[18:40].strip() else None,
             "id": dbref[47:67].strip() if dbref[42:54].strip() else None,
             "db_seq_begin": int(dbref2[45:55].strip()) if dbref2[45:55].strip() else None,
             "db_seq_end": int(dbref[57:67].strip()) if dbref2[57:67].strip() else None
            }

        #Process SEQADV
        seqadv = self.get_records_by_name("SEQADV")
        self.sequence_differences = [{
         "code": s[7:11] if s[7:11].strip() else None,
         "res_name": s[12:15].strip() if s[12:15].strip() else None,
         "chain": s[16] if s[16].strip() else None,
         "seq_num": int(s[18:22].strip()) if s[18:22].strip() else None,
         "i_code": s[22] if s[22].strip() else None,
         "database": s[24:28].strip() if s[24:28].strip() else None,
         "accession": s[30:37].strip() if s[30:37].strip() else None,
         "db_res": s[39:42].strip() if s[39:42].strip() else None,
         "db_seq": int(s[43:48].strip()) if s[43:48].strip() else None,
         "comment": s[49:70].strip() if s[49:70].strip() else None
        } for s in seqadv]

        #Process SEQRES
        seqres = self.get_records_by_name("SEQRES")
        chains = list(set([r[11] for r in seqres]))
        chains.sort()
        self.sequences = []
        for chain in chains:
            residues = " ".join([r[19:].strip() for r in seqres if r[11] == chain]
             ).replace("  ", " ").split()
            self.sequences.append({
             "chain": chain,
             "residues": residues
            })

        #Process MODRES
        modres = self.get_records_by_name("MODRES")
        self.residue_modifications = [{
         "code": r[7:11] if r[7:11].strip() else None,
         "res_name": r[12:15].strip() if r[12:15].strip() else None,
         "chain_id": r[16] if r[16].strip() else None,
         "seq_num": int(r[18:22].strip()) if r[18:22].strip() else None,
         "i_code": r[22] if r[22].strip() else None,
         "std_res": r[24:27].strip() if r[24:27].strip() else None,
         "comment": r[29:70].strip() if r[29:70].strip() else None
        } for r in modres]



class HeterogenSection(PdbSection):

    RECORD_NAMES = ("HET", "HETNAM", "HETSYN", "FORMUL")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process HETs
        hets = self.get_records_by_name("HET")
        self.hets = [{
         "het_id": h[7:10].strip() if h[7:10] else None,
         "chain": h[12] if h[12].strip() else None,
         "seq_num": int(h[13:17].strip()) if h[13:17].strip() else None,
         "insert": h[17] if h[17].strip() else None,
         "num_atoms": int(h[20:25].strip()) if h[20:25].strip() else None,
         "description": h[30:70].strip() if h[30:70].strip() else None
        } for h in hets]

        #Process HETNAMs
        hetnams = self.get_records_by_name("HETNAM")
        names = set([h[11:14] for h in hetnams])
        self.hetnams = []
        for name in names:
            fullname = " ".join([h[15:].strip() for h in hetnams if h[11:14] == name]
             ).replace("  ", " ")
            self.hetnams.append({
             "code": name,
             "fullname": fullname
            })

        #Process HETSYNs
        hetsyns = self.get_records_by_name("HETSYN")
        names = set([h[11:14] for h in hetsyns])
        self.hetsyns = []
        for name in names:
            synonyms = " ".join([h[15:].strip() for h in hetsyns if h[11:14] == name]
             ).replace("  ", " ").replace(", ", ",").split(",")
            self.hetsyns.append({
             "code": name,
             "synonyms": synonyms
            })

        #Process FORMULs
        formuls = self.get_records_by_name("FORMUL")
        names = set([h[12:15].strip() for h in formuls])
        self.formuls = []
        for name in names:
            records = [r for r in formuls if r[12:15].strip() == name]
            self.formuls.append({
             "component_number": int(records[0][8:10].strip()) if records[0][8:10].strip() else None,
             "het_id": name,
             "water": records[0][18] == "*",
             "formula": " ".join([r[19:].strip() for r in records if r[12:15].strip() == name]
              ).replace("  ", " ")
            })



class SecondaryStructureSection(PdbSection):

    RECORD_NAMES = ("HELIX", "SHEET")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process HELIXs
        helices = self.get_records_by_name("HELIX")
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
        sheet_lines = self.get_records_by_name("SHEET")
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



class ConnectivityAnnotationSection(PdbSection):

    RECORD_NAMES = ("SSBOND", "LINK", "CISPEP")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process SSBONDs
        ssbonds = self.get_records_by_name("SSBOND")
        self.ssbonds = [{
         "serial_num": int(s[7:10].strip()) if s[7:10].strip() else None,
         "residue_1_name": s[11:14].strip() if s[11:14].strip() else None,
         "residue_1_chain": s[15] if s[15].strip() else None,
         "residue_1_number": int(s[17:21].strip()) if s[17:21].strip() else None,
         "residue_1_insert": s[21] if s[21].strip() else None,
         "residue_1_symmetry": s[59:65].strip() if s[59:65].strip() else None,
         "residue_2_name": s[25:28].strip() if s[25:28].strip() else None,
         "residue_2_chain": s[29] if s[29].strip() else None,
         "residue_2_number": int(s[31:35].strip()) if s[31:35].strip() else None,
         "residue_2_symmetry": s[66:72].strip() if s[59:65].strip() else None,
         "residue_2_insert": s[35] if s[35].strip() else None,
         "disulfide_distance": float(s[73:78].strip()) if s[73:78].strip() else None
        } for s in ssbonds]

        #Process LINKs
        links = self.get_records_by_name("LINK")
        self.links = [{
         "residue_1_atom": s[12:16].strip() if s[12:16].strip() else None,
         "residue_1_name": s[17:20].strip() if s[17:20].strip() else None,
         "residue_1_chain": s[21] if s[21].strip() else None,
         "residue_1_number": int(s[22:26].strip()) if s[22:26].strip() else None,
         "residue_1_symmetry": s[59:65].strip() if s[59:65].strip() else None,
         "residue_1_insert": s[21] if s[21].strip() else None,
         "residue_2_atom": s[42:46].strip() if s[42:46].strip() else None,
         "residue_2_name": s[47:50].strip() if s[47:50].strip() else None,
         "residue_2_chain": s[51] if s[51].strip() else None,
         "residue_2_number": int(s[52:56].strip()) if s[52:56].strip() else None,
         "residue_2_symmetry": s[66:72].strip() if s[59:65].strip() else None,
         "link_distance": float(s[73:78].strip()) if s[73:78].strip() else None
        } for s in links]

        #Process CISPEPs
        cispeps = self.get_records_by_name("CISPEP")
        self.cispeps = [{
         "serial_num": int(s[7:10].strip()) if s[7:10].strip() else None,
         "residue_1_name": s[11:14].strip() if s[11:14].strip() else None,
         "residue_1_chain": s[15] if s[15].strip() else None,
         "residue_1_number": int(s[17:21].strip()) if s[17:21].strip() else None,
         "residue_1_insert": s[21] if s[21].strip() else None,
         "residue_2_name": s[25:28].strip() if s[25:28].strip() else None,
         "residue_2_chain": s[29] if s[29].strip() else None,
         "residue_2_number": int(s[31:35].strip()) if s[31:35].strip() else None,
         "residue_2_insert": s[35] if s[35].strip() else None,
         "modnum": int(s[43:46].strip()) if s[43:46].strip() else None,
         "angle_measure": float(s[53:59].strip()) if s[53:59].strip() else None
        } for s in cispeps]



class MiscellaneousSection(PdbSection):

    RECORD_NAMES = ("SITE")

    def __init__(self, *args, **kwargs):
        PdbSection.__init__(self, *args, **kwargs)

        #Process SITEs
        sites = self.get_records_by_name("SITE")
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
        cryst1 = self.get_records_by_name("CRYST1")
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
        origx1 = self.get_records_by_name("ORIGX1")
        origx2 = self.get_records_by_name("ORIGX2")
        origx3 = self.get_records_by_name("ORIGX3")
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
        scale1 = self.get_records_by_name("SCALE1")
        scale2 = self.get_records_by_name("SCALE2")
        scale3 = self.get_records_by_name("SCALE3")
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
        mtrix1 = self.get_records_by_name("MTRIX1")
        mtrix2 = self.get_records_by_name("MTRIX2")
        mtrix3 = self.get_records_by_name("MTRIX3")
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
