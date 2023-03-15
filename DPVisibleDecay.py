from DPInvisibleDecay import DPInvisibleDecay

class DPVisibleDecay(DPInvisibleDecay):
    'Calculate yield of visible decay of dark photon, draw the exclusion region'

    def __init__(self) -> None:
        super().__init__()

        self.EOT = 1e14
        self.E0 = 8