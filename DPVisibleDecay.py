from DPInvisibleDecay import DPInvisibleDecay
from scipy.integrate import dblquad
import numpy as np

class DPVisibleDecay(DPInvisibleDecay):
    'Calculate yield of visible decay of dark photon, draw the exclusion region'

    def __init__(self, EOT, E0, thickness_in_X0=0.1, verbose=0):
        super().__init__(EOT, E0, thickness_in_X0, verbose)

        self.signal_efficiency = 0.4
        self.bkg = 10
        self.detect_range = [1, 17] # in cm

    def SetSignalEfficiency(self, signal_efficiency):
        self.signal_efficiency = signal_efficiency

    # Set detect range, 0 point is front end of target (in cm)
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
        thickness_in_cm = self.thickness / self.rho
        EA = x * self.E0
        z0, z1 = self.detect_range
        length = self.GetDecayLength(EA, mA, epsilon)
        
        # Consider that target has a thickness anymore, in case of small decay length,
        # the position of decay matters
        if (length > 10 * thickness_in_cm): 
            thickness_effect = 1
        elif (length < 0.1 * thickness_in_cm):
            thickness_effect = 0
        else:
            fraction = thickness_in_cm / length
            thickness_effect = 1 / fraction * (np.exp(fraction) - 1)

        dXS_dx = super().GetdXS_dx(x, mA, epsilon)
        return dXS_dx * (np.exp(- z0 / length) - np.exp(- z1 / length)) * thickness_effect

    # Get acceptance of signal at each point
    def GetSignalEfficiency(self, mA, epsilon):
        return self.signal_efficiency

    # For cross check with paper
    # def GetDecayLength(self, EA, mA, epsilon):
    #     alpha = self.alpha
    #     me = self.me
    #     mmu = self.mmu
    #     # Only consider A'->ee and A'->mumu channel
    #     Neff = 0
    #     if (mA > 2 * me):
    #         Neff = 1
    #     if (mA > 2 * mmu):
    #         Neff = 2
    #     length =  80 / Neff * (1e-4 / epsilon)**2 * (0.1 / mA) - 349 / EA * mA
    #     if (length > 1):
    #         return length
    #     else:
    #         return 1
    