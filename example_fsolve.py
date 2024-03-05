from scipy.optimize import fsolve
from DPVisibleDecay import DPVisibleDecay
import matplotlib.pyplot as plt
import numpy as np

# Initialize and configure ExclusionDrawer
darkSHINE1 = DPVisibleDecay(3e14, 8, verbose=0)
darkSHINE1.SetBkgYield(10)
darkSHINE1.initialization()

efficiency = 0.6
darkSHINE1.SetSignalEfficiency(efficiency)

minrange = 5
darkSHINE1.SetDetectRange([minrange, 17])

# Define the mA points to calculate epsilon.
# This need to be adjusted mannually for more scan points in sensitive region.
mAs_1 = 10 ** np.arange(-2.8, -1.4, 0.1)
mAs_2 = 10 ** np.arange(-1.4, -0.9, 0.01)
mAs = np.append(mAs_1, mAs_2)

# Define the equation to slove for 90% CL
def equation(x):
    return darkSHINE1.GetSignificance(mA, x) - 1.64

# Store the result
xs, ys_min, ys_max = [], [], []

# Here begin to solve the equation numerically.
# Start with large root to get upper edge of exclusion limit
init_root = 0.01
for count, mA in enumerate(mAs):
    result = fsolve(equation, x0=init_root, full_output=True)
    root, status = result[0][0], result[2]
    scale = 1 if count == 0 else root / ys_max[count - 1]
    print(f" ({count+1}/{len(mAs)}) {root = :.4e}, {status = }, {scale  =}")
    if status == 1:
        xs.append(mA)
        ys_max.append(root)

        # The initial root is given by last scan point
        # If the difference between scan point is large, a scale for initial root is given manually.
        init_scale = 0.75 if count < len(mAs_1) else 1
        init_root = root * init_scale 
    else:
        break

# Start with small root to get lower edge of exclusion limit
init_root = 0.0001
for count, mA in enumerate(mAs):
    result = fsolve(equation, x0=init_root, full_output=True)
    root, status = result[0][0], result[2]
    print(f" ({count+1}/{len(mAs)}) {root = :.4e}, {status = }")
    if status == 1:
        ys_min.append(root)
        init_root = root
    else:
        break

# Store the data
with open(f"DarkShine_minRange5_eff{int(efficiency*100)}.lmt", 'w') as f:
    f.write("#      mass       lower       upper\n")
    for i in range(len(xs)):
        f.write(f" {xs[i]:.4e}  {ys_min[i]:.4e}  {ys_max[i]:.4e}\n")
    f.close()