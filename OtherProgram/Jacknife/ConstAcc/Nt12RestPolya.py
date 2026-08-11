import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads = "H:\\ConstAccLarge\\Nt12\\Polyakov\\"
folderHeads2 = "H:\\ConstAccLarge\\Nt12\\Polyakov2\\"

lst = ["68", "69", "70", "71", "72", "73", "74", "75", "76"]

p = []
pe = []
ps = []
pse = []
pt = []

for i in range(len(lst)):
    fileNames = folderHeads + "BSQ__{}_polyakov.csv".format(lst[i])
    testarray = LoadMathematicaCSV(fileNames)
    print(np.shape(testarray))
    fileNames = folderHeads2 + "BSQ__{}_polyakov.csv".format(lst[i])
    testarray = np.append(testarray, LoadMathematicaCSV(fileNames))
    print(np.shape(testarray))
    testarrayre = np.abs(testarray)
    v1, s1 = JacknifeMean(testarrayre, lst[i])
    v2, s2 = JacknifeCumulant(testarrayre, lst[i])
    _, _, t = AutoCorrelationSingleVariable(testarrayre)
    p.append(v1)
    pe.append(s1)
    ps.append(v2)
    pse.append(s2)
    pt.append(t)


p = np.array(p) / 3
pe = np.array(pe) / 3
ps = 12 * 12 * 12 * 24 * 24 * 7 * np.array(ps) / 9
pse = 12 * 12 * 12 * 24 * 24 * 7 * np.array(pse) / 9
pt = np.array(pt)
pe = pe * np.sqrt(2 * pt)
pse = pse * np.sqrt(2 * pt)
print(PrintAsMatlabArray(p, "polya"))
print(PrintAsMatlabArray(pe, "polyae"))
print(PrintAsMatlabArray(ps, "polyaSusp"))
print(PrintAsMatlabArray(pse, "polyaeSusp"))
print(PrintAsMatlabArray(pt, "polyat"))

xlst = [278.1, 289.1, 303.0, 307.0, 310.9, 323.9, 333.3, 339.3, 349.8]
errorbar(xlst, [ps], [pse], xlabel="$T({\\rm MeV})$", ylabel="$\chi _L/T^3$")

