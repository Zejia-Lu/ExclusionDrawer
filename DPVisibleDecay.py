from DPInvisibleDecay import DPInvisibleDecay
from scipy.integrate import dblquad
import numpy as np

class DPVisibleDecay(DPInvisibleDecay):
    'Calculate yield of visible decay of dark photon, draw the exclusion region'

    def __init__(self, EOT, Ee, verbose=0):
        super().__init__(EOT, Ee, verbose)

        self.bkg = 10
        self.detect_range = [1, 18] # in cm

    # Set detect range (in cm)
    def SetDetectRange(self, detect_range):
        self.detect_range = detect_range

    # Get decay length in cm
    def GetDecayLength(self, EA, mA, epsilon):
        alpha = self.alpha
        me = self.me
        mmu = self.mmu
        # Only consider A'->ee and A'->mumu channel
        Gamma_ee, Gamma_mumu = 0., 0.
        if (mA > 2 * me):
            Gamma_ee = 1/3 * alpha * epsilon**2 * mA * (1 - 4 * me**2 / mA**2)**0.5 * (1 + 2 * me**2 / mA**2)
        if (mA > 2 * mmu):
            Gamma_mumu = 1/3 * alpha * epsilon**2 * mA * (1 - 4 * mmu**2 / mA**2)**0.5 * (1 + 2 * mmu**2 / mA**2)
        length =  1 / (Gamma_ee + Gamma_mumu) * EA / mA / self.cm
        return length
    
    # Get differential cross section, has integrated on z
    def GetdXS_dx(self, x, mA, epsilon):
        EA = x * self.E0
        z0, z1 = self.detect_range
        length = self.GetDecayLength(EA, mA, epsilon)
        dXS_dx = super().GetdXS_dx(x, mA, epsilon)
        return dXS_dx / length**2 * (np.exp(- z0 / length) - np.exp(- z1 / length))

    # Get acceptance of signal at each point
    def GetSignalEfficiency(self, mA, epsilon):
        return 0.4
