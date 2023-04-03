from DPVisibleDecay import DPVisibleDecay
import matplotlib.pyplot as plt
import numpy as np

darkSHINE1 = DPVisibleDecay(3e14, 8, verbose=1)
darkSHINE1.SetBkgYield(10)

# Set mA and epsilon to be calculated (option)
mAs = np.array([2., 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500], dtype='double')
mAs *= 1e-3
epsilons = np.array([1e-5, 2e-5, 3e-5, 5e-5, 8e-5, 1e-4, 2e-4, 3e-4, 6e-4, 9e-4, 1e-3, 2e-3, 4e-3, 6e-3, 1e-2])
darkSHINE1.SetmAs(mAs)
darkSHINE1.SetEpsilons(epsilons)

# Initialize to avoid repeated calculation
darkSHINE1.initialization()

# Calculate grid and draw contour plot
X,Y,Z = darkSHINE1.GetGrid()

plt.contour(X, Y, Z, levels=[0.9])
plt.xlabel('$log_{10}(m_A)$')
plt.ylabel('$log_{10}\epsilon$')
plt.show()