from ExclusionDrawer import ExclusionDrawer
from scipy.integrate import quad
from math import log
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

class DPInvisibleDecay(ExclusionDrawer):
    'Calculate yield of invisible decay of dark photon, draw the exclusion region'

    def __init__(self, EOT, E0, thickness_in_X0=0.1, verbose=0):
        super().__init__(EOT, E0, thickness_in_X0, verbose)

        self.bkg = 0.015

        self.mAs = np.array([1., 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 700, 1000], dtype='double') * 1e-3
        self.epsilons = 10 ** np.arange(-7, -1, 0.5)
    
    # Set background yield
    def SetBkgYield(self, bkg):
        self.bkg = bkg

    def SetmAs(self, mAs):
        self.mAs = mAs
    
    def SetEpsilons(self, epsilons):
        self.epsilons = epsilons

    # Initialize to avoid repeated calculation
    def initialization(self):
        Z = DPInvisibleDecay.Z
        A = DPInvisibleDecay.A
        me = DPInvisibleDecay.me

        self.a_el = 111 * Z**(-1/3) / me
        self.a_in = 773 * Z**(-2/3) / me
        self.d = 0.164 * A**(-2/3)

    # Calculate the effective flux of photons according to (A17) in [arXiv:0906.0580]
    def GetChi(self, mA):
        Z = DPInvisibleDecay.Z
        mp = DPInvisibleDecay.mp
        mup = DPInvisibleDecay.mup
        E0 = self.E0
        a_el = self.a_el
        a_in = self.a_in
        d = self.d

        def G2_el(t):
            return (a_el**2 * t / (1 + a_el**2 * t))**2 * (1 / (1 + t / d))**2 * Z**2
        
        def G2_in(t):
            return (a_in**2 * t / (1 + a_in**2 * t))**2 * ((1 + t / 4 / mp**2 * (mup**2 - 1)) / (1 + t / 0.71)**4)**2 * Z
        
        t_min = (mA**2 / 2 / E0)**2
        t_max = mA**2

        def dChi_dt(t):
            return (t - t_min) / t**2 * (G2_el(t) + G2_in(t))
            # return (t - t_min) / t**2 * G2_el(t)
        
        Chi, error = quad(dChi_dt, t_min, t_max)
        return Chi
    
    def GetdChi_dt(self, t, mA):
        Z = DPInvisibleDecay.Z
        mp = DPInvisibleDecay.mp
        mup = DPInvisibleDecay.mup
        E0 = self.E0
        a_el = self.a_el
        a_in = self.a_in
        d = self.d

        def G2_el(t):
            return (a_el**2 * t / (1 + a_el**2 * t))**2 * (1 / (1 + t / d))**2 * Z**2
        
        def G2_in(t):
            return (a_in**2 * t / (1 + a_in**2 * t))**2 * ((1 + t / 4 / mp**2 * (mup**2 - 1)) / (1 + t / 0.71)**4)**2 * Z
        
        t_min = (mA**2 / 2 / E0)**2
        t_max = mA**2

        def dChi_dt(t):
            return (t - t_min) / t**2 * (G2_el(t) + G2_in(t))
            # return (t - t_min) / t**2 * G2_el(t)
        
        return dChi_dt(t)

    # Calculate differential cross section in pb according to (A14) in [arXiv:0906.0580]
    def GetdXS_dx(self, x, mA, epsilon):
        alpha = DPInvisibleDecay.alpha
        me = DPInvisibleDecay.me
        E0 = self.E0

        return 4 * alpha**3 * epsilon**2 * self.GetChi(mA) * (1 - mA**2 / E0**2)**0.5 * (1 - x + 1/3 * x**2) / (mA**2 * (1 - x) / x + me**2 * x) * self.toPB

    # Calculate total cross section by integrate over differential cross section
    def GetXS(self, mA, epsilon):
        def GetdXS_dx_at_point(x):
            return self.GetdXS_dx(x, mA, epsilon)
        XS, error = quad(GetdXS_dx_at_point, mA/self.E0, 1 - mA/self.E0)
        # print(f"{mA = }, {XS = }, {error = }")
        return XS
    
    # Calculate total cross section by formula according to (A15) in [arXiv:0906.0580]
    def GetXS_formula(self, mA, epsilon):
        alpha = DPInvisibleDecay.alpha
        me = DPInvisibleDecay.me
        E0 = self.E0

        beta = (1 - mA**2 / E0**2)**0.5
        critical = max(me**2 / mA**2, mA**2 / E0**2)

        return 4. / 3. * alpha**3 * epsilon**2 * beta / mA**2 * self.GetChi(mA) * np.log10(1 / critical) * self.toPB

    # Get acceptance of signal at each point
    def GetSignalEfficiency(self, mA, epsilon):
        return 0.5

    # Calculate signal yield from XS and beam/target information
    def GetSignalYield(self, mA, epsilon):
        return self.GetSignalEfficiency(mA, epsilon) * self.GetXS(mA, epsilon) * self.thickness * self.EOT * self.NA / self.A * 1e-36 
    
    # Get significance Z0 from Asimov's formula
    def GetSignificance(self, mA, epsilon):
        s = self.GetSignalYield(mA, epsilon)
        b = self.bkg
        Z2 = 2 * ((s + b) * log(1 + s / b) - s)
        if (Z2 > 0):
            Z = Z2 ** 0.5
        else :
            Z = 0.
        
        if (self.verbose > 0):
            print('{:0} {:1} {:2}'.format(mA, epsilon, Z))
        
        return Z

    # Calculate signal yield for each point, just for cross check
    def GetYieldGrid(self):
        X, Y = np.meshgrid(self.mAs, self.epsilons)
        Z = X + Y
        for i in range(len(X)):
            for j in range(len(X[0])):
                Z[i][j] = self.GetSignalYield(X[i][j], Y[i][j])
                if (self.verbose > 0):
                    print('{:0} {:1} {:2}'.format(X[i][j], Y[i][j], Z[i][j]))

        X = np.log10(X)
        Y = np.log10(Y)
        return X, Y, Z

    # Calculate confidence level for each point
    def GetGrid(self):
        X, Y = np.meshgrid(self.mAs, self.epsilons)
        Z = X + Y
        for i in range(len(X)):
            for j in range(len(X[0])):
                Z[i][j] = norm.cdf(self.GetSignificance(X[i][j], Y[i][j]))

        X = np.log10(X)
        Y = np.log10(Y)
        return X, Y, Z

    # Draw exclusion plot
    def DrawExclusion(self, confidence):
        X, Y, Z = self.GetGrid()
        plt.rc('text', usetex=True)
        plt.contour(X, Y, Z, levels=[confidence])
        plt.xlabel('$log_{10}(m_A)$')
        plt.ylabel('$log_{10}\epsilon$')
        plt.show()
