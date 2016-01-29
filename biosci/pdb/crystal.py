class UnitCell:

    def __init__(self, crystal_section):
        self.a = crystal_section.a
        self.b = crystal_section.b
        self.c = crystal_section.c
        self.alpha = crystal_section.alpha
        self.beta = crystal_section.beta
        self.gamma = crystal_section.gamma
        self.s_group = crystal_section.s_group
        self.z = crystal_section.z


    def __repr__(self):
        return "<Unit Cell>%i×%i×%i" % (self.a, self.b, self.c)
