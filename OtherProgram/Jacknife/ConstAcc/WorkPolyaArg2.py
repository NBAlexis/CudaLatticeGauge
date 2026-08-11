import numpy as np
from matplotlib import pyplot as plt

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMatlabArray

lst = ["00", "05", "10", "15", "20", "25", "30", "35", "40"]
gv = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04]

args = []
argse = []
c4s = []
c4se = []
for i in range(len(lst)):
    f1 = "G:\\ConstAccNew\\60\\Polyakov\\ACCQ60__{}_polyakov.csv".format(lst[i])
    f2 = "G:\\ConstAccNew\\60\\Chiral\\ACCQ60__{}_condensateZSlicepCCLightCMTKSGamma4.csv".format(lst[i])
    arr1 = np.angle(LoadMathematicaCSV(f1))
    arr2 = -np.imag(LoadMathematicaCSV(f2))[:, 0] / 2
    v, s = JacknifeMean(arr1)
    _, _, t = AutoCorrelationSingleVariable(arr1)
    args.append(v)
    argse.append(s * np.sqrt(2 * t))
    v, s = JacknifeMean(arr2)
    _, _, t = AutoCorrelationSingleVariable(arr2)
    c4s.append(v)
    c4se.append(s * np.sqrt(2 * t))

print(PrintAsMatlabArray(args, "args"))
print(PrintAsMatlabArray(argse, "argse"))
print(PrintAsMatlabArray(c4s, "c4s"))
print(PrintAsMatlabArray(c4se, "c4se"))

plt.plot(gv, args)
plt.show()
plt.plot(gv, c4s)
plt.show()