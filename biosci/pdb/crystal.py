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




class Transformation:

    def __init__(self, crystal_section, m, n):
        for a in range(1,4):
            for b in range(1,4):
                self.__dict__["%s%i%i" % (m, a, b)] = crystal_section.__dict__["%s%i%i" % (m, a, b)]
            self.__dict__["%s%i" % (n, a)] = crystal_section.__dict__["%s%i" % (n, a)]



class SubmittedCoordinatesTransformation(Transformation):

    m, n = "o", "t"

    def __init__(self, crystal_section, *args, **kwargs):
        Transformation.__init__(self, crystal_section, self.m, self.n, *args, **kwargs)



class CrystallographicCoordinatesTransformation(Transformation):

    m, n = "s", "u"

    def __init__(self, crystal_section, *args, **kwargs):
        Transformation.__init__(self, crystal_section, self.m, self.n, *args, **kwargs)



class MatrixTransformation(Transformation):

    m, n = "m", "v"

    def __init__(self, crystal_section, *args, **kwargs):
        Transformation.__init__(self, crystal_section, self.m, self.n, *args, **kwargs)
        self.serial = crystal_section.serial
        self.i_given = crystal_section.i_given
