import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray, \
    PrintAsMatlabArray
from Visualization import errorbar

paramlst = ["0521", "10209", "2020"]
kappalst = [0.21, 0.209, 0.2]
header = "I:\\StrongCoupling\\WilsonDirac\\Chiral\\"

def getPrAndNumberOfSites():
    Lx = 14
    nEdge = (Lx - 2) // 2
    numSite1 = np.array([0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    pr1 = np.array([0.0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    for x in range(-nEdge, nEdge + 1):
        for y in range(-nEdge, nEdge + 1):
            prSq = x * x + y * y
            pr1[prSq] = prSq
            numSite1[prSq] = numSite1[prSq] + 1
    pr2 = pr1[numSite1 > 0]
    numSite2 = numSite1[numSite1 > 0]
    idxedge = np.where(pr2 == nEdge * nEdge)[0][0]
    pr2 = np.sqrt(pr2)
    return pr2, numSite2, np.sum(numSite2), idxedge

def getMean(filename: str, kappa: float):
    resvAll = []
    ressAll = []
    resvIn = []
    ressIn = []
    _, nums, sumnum, edgenum = getPrAndNumberOfSites()
    for i in range(0, 11):
        # NC0521__condensateChiral_Nt3_O0
        arr = LoadMathematicaCSV(f"{header}NC{filename}__condensateChiral_Nt3_O{i}.csv")
        rearAll = np.dot(np.real(arr), nums) / sumnum
        v, s = JacknifeMean(rearAll)
        _, _, t = AutoCorrelationSingleVariable(rearAll)
        resvAll.append(v)
        ressAll.append(s * np.sqrt(2 * t))
        rearIn = np.dot(np.real(arr)[:,:edgenum], nums[:edgenum]) / np.sum(nums[:edgenum])
        v, s = JacknifeMean(rearIn)
        _, _, t = AutoCorrelationSingleVariable(rearIn)
        resvIn.append(v)
        ressIn.append(s * np.sqrt(2 * t))
    factor = 2 * kappa
    return np.vstack((factor * np.array(resvIn),
                      factor * np.array(ressIn),
                      factor * np.array(resvAll),
                      factor * np.array(ressAll)))

finalres = None
for paramidx in range(len(paramlst)):
    if finalres is None:
        finalres = getMean(paramlst[paramidx], kappalst[paramidx])
    else:
        finalres = np.vstack((finalres, getMean(paramlst[paramidx], kappalst[paramidx])))

np.savetxt("Data/ChiralWD.csv", finalres, delimiter=",")
print(finalres)



