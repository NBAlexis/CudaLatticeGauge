import numpy as np
from matplotlib import pyplot as plt

from JacknifePrograms import LoadMathematicaCSV


# f1 = "G:\\ConstAccNewMid\\595\\Polyakov\\ACCQM595__10_polyakov.csv"
# f2 = "G:\\ConstAccNewMid\\60\\Polyakov\\ACCQM60__10_polyakov.csv"
# f3 = "G:\\ConstAccNewMid\\605\\Polyakov\\ACCQM605__10_polyakov.csv"
f1 = "G:\\ConstAccNewMid\\595\\Polyakov\\ACCQM595__06_polyakov.csv"
f2 = "G:\\ConstAccNewMid\\595\\Polyakov\\ACCQM595__07_polyakov.csv"
f3 = "G:\\ConstAccNewMid\\605\\Polyakov\\ACCQM605__00_polyakov.csv"

print(np.mean(np.abs(LoadMathematicaCSV(f1))))
print(np.mean(np.abs(LoadMathematicaCSV(f3))))
arr1 = np.angle(LoadMathematicaCSV(f1))
arr2 = np.angle(LoadMathematicaCSV(f2))
arr3 = np.angle(LoadMathematicaCSV(f3))

print(np.mean(arr1), np.mean(arr2), np.mean(arr3))
# print(np.std(arr1), np.std(arr2), np.std(arr3))

fig = plt.figure(figsize=(10, 6))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal Axes and the main Axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(1, 2,  width_ratios=(4, 1),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
# Create the Axes.
ax = fig.add_subplot(gs[0, 0])
ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)

xpoints = [i + 101 for i in range(9900)]
ax.scatter(xpoints, arr1, s=1, label='ag=0')
ax.scatter(xpoints, arr2, s=1, label='ag=0.06')
ax.scatter(xpoints, arr3, s=1, label='ag=0.07')
ax.set_xlabel("TU")
ax.set_ylabel("Im($c_4(z)$)")
ax.legend()

ax_histy.hist(arr1, bins=50, range=[-0.05, 0.05], histtype='step', orientation='horizontal')
ax_histy.hist(arr2, bins=50, range=[-0.05, 0.05], histtype='step', orientation='horizontal')
ax_histy.hist(arr3, bins=50, range=[-0.05, 0.05], histtype='step', orientation='horizontal')
plt.show()
# plt.savefig('c4scatter2.eps')

# print(np.mean(arr2), np.mean(arr2[:4000]), np.mean(arr2[2000:]))
