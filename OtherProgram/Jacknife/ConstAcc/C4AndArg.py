import numpy as np
from matplotlib import pyplot as plt

from JacknifePrograms import LoadMathematicaCSV

f1 = "H:\\ConstAccNewMid\\585\\Chiral\\ACCQM585__07_condensateZSlicepCCLightCMTKSGamma4.csv"
f2 = "H:\\ConstAccNewMid\\585\\Polyakov\\ACCQM585__07_polyakov.csv"

fig = plt.figure(figsize=(10, 6))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal Axes and the main Axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 1, height_ratios=(1, 1),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_arg = fig.add_subplot(gs[1, 0], sharex=ax)

xpoints = [i + 101 for i in range(9900)]
arr1 = -np.imag(LoadMathematicaCSV(f1)[:, 0]) / 2
arr2 = np.angle(LoadMathematicaCSV(f2))
ax.scatter(xpoints, arr1, s=1)
ax.set_ylabel("-Im($c_4(z)$)")

ax_arg.scatter(xpoints, arr2, s=1)
ax_arg.set_ylabel("Arg(L)")
ax_arg.set_xlabel("TU")

plt.savefig('c4arg2.eps')
