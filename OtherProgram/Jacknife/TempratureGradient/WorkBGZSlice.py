import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray

foldername = "G:/BetaGradient/BetaScan/BS529536/Chiral/BS02__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:/BetaGradient/BetaGradient/{}/Polyakov{}/{}__{}_polyakov_ZSlice.csv"
# foldername = "G:/BetaGradient/BetaGradient/{}/Chiral{}/{}__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:/BetaGradient/BetaGradient/{}/Chiral{}/{}__{}_condensateZSlicepCCLightCMTKSGamma4.csv"
# foldername = "G:/BetaGradient/BetaGradient/{}/Chiral{}/{}__{}_condensateZSlicepCCLightChiralKS.csv"
# lst = ["529", "530", "531", "532", "533", "534", "535", "536"]
lst = ["529", "536"]
# lst2 = ["BG005", "BG01", "BG015", "BG02"]
# lst2 = ["BG005", "BG02"]
lst2 = [""]
# lst2 = [""]
# lst3 = [1, 2]

# signs = [[1, 1], [1, 1], [1, 1], [1, 1]]
# signs = [[1, -1], [-1, -1], [1, -1], [1, -1]]
signs = [[1, -1], [1, -1]]

for i in range(len(lst2)):
    for j in range(len(lst)):
        polyalst = []
        polyalste = []
        susplst = []
        susplste = []
        times = []
        runname = "{}-{}".format(lst2[i], lst[j])
        filename = foldername.format(lst[j])
        # filename = foldername.format(lst[j])
        testarray = LoadMathematicaCSV(filename)
        # filename = foldername.format(lst2[i], lst3[1], lst2[i], lst[j])
        # testarray = np.vstack((testarray * signs[i][0], LoadMathematicaCSV(filename) * signs[i][1]))
        # arr = np.real((testarray[:, 0] + testarray[:, 6])/2)

        for k in range(0, 12):
            arr = -np.real(testarray[:, k])
            v, s = JacknifeMean(arr)
            _, _, t = AutoCorrelationSingleVariable(arr)
            # print(v)
            # print(s)
            polyalst.append(v)
            polyalste.append(s)
        # v, s = JacknifeCumulant(arr, "susp" + runname)
        # print(v)
        # print(s)
        # susplst.append(v)
        # susplste.append(s)
            times.append(t)
    # print(PrintAsMathematicaArray(polyalst, "chiral" + lst2[i]))
    # print(PrintAsMathematicaArray(polyalste, "chirale" + lst2[i]))
    # print(PrintAsMathematicaArray(susplst, "chiralsusp" + lst2[i]))
    # print(PrintAsMathematicaArray(susplste, "chiralsuspe" + lst2[i]))
    # print(PrintAsMathematicaArray(times, "time" + lst2[i]))
    # print(PrintAsMathematicaArray(polyalst, "polya" + lst2[i]))
    # print(PrintAsMathematicaArray(polyalste, "polyae" + lst2[i]))
    # print(PrintAsMathematicaArray(susplst, "polyasusp" + lst2[i]))
    # print(PrintAsMathematicaArray(susplste, "polyasuspe" + lst2[i]))
    # print(PrintAsMathematicaArray(times, "ptime" + lst2[i]))
        print(PrintAsMathematicaArray(polyalst, "polya" + lst2[i] + lst[j]))
        print(PrintAsMathematicaArray(polyalste, "polyae" + lst2[i] + lst[j]))
        print(PrintAsMathematicaArray(times, "ptime" + lst2[i] + lst[j]))
