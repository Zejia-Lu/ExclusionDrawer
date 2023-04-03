class ExclusionDrawer:
    'Draw exclusion region, to be inherited from specified model'

    # Properties of W target
    Z = 74
    A = 184
    X0 = 6.76 # in g/cm^2
    rho = 19.35 # in g/cm^3
    
    # Constant
    alpha = 1. / 137
    me = 0.511e-3 # in GeV
    mmu = 0.10566 # in GeV
    NA = 6.02e23
    mp = 0.938
    mup = 2.79
    cm = 5.0677e13 # in GeV^-1
    toPB = 1e36 / cm**2 # from GeV^-2 to pb

    def __init__(self, EOT, E0, thickness_in_X0=0.1, verbose=0):
        self.EOT = EOT
        self.E0 = E0
        self.verbose = verbose
        self.thickness = self.X0 * thickness_in_X0

    # Set method
    def SetEOT(self, EOT):
        self.EOT = EOT
    
    def SetE0(self, E0):
        self.E0 = E0
    
    def SetVerbose(self, verbose):
        self.verbose = verbose

    def GetXS(self):
        pass

    def GetSignalYield(self):
        pass

    def GetBkgYield(self):
        pass

    def GetSignificance(self):
        pass

    def DrawExclusion(self):
        pass