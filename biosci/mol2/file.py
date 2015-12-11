class Mol2File:
    """A representation of the mol2 file itself, not the structure it represents."""

    def __init__(self, mol2_contents):
        lines = [line for line in mol2_contents.split("\n")]

        self.records = []
        for line_number, line in enumerate(lines, start=1):
            if not line.strip():
                self.records.append(BlankRecord(line, line_number))
            elif line.strip()[0] == "#":
                self.records.append(CommentRecord(line, line_number))
            elif line.strip()[0] == "@":
                self.records.append(RtiRecord(line, line_number))
            else:
                self.records.append(DataRecord(line, line_number))



class Record:
    """A mol2 record (a line in the file)."""

    def __init__(self, line, line_number):
        self.line_number = line_number
        self.text = line


    def __repr__(self):
        return "<Record %i>%s" % (self.line_number, self.text)


class BlankRecord(Record):
    pass



class CommentRecord(Record):
    pass



class RtiRecord(Record):
    pass



class DataRecord(Record):
    pass
