import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray, PrintAsMatlabArray

# foldername = "H:\\BGQuench\\Polyakov\\BGQ__{}_polyakov_ZSlice.csv"
# foldername = "H:\\BGQuench\\Chiral\\BGQ__{}_condensateZSlicepCCLightChiralKS.csv"
# lst = ["570", "575", "580", "585", "590", "595", "600", "605"]

# foldername = "H:\\BetaGradient\\BetaScan\\BS529536\\Polyakov\\BS02__{}_polyakov_ZSlice.csv"
# foldername = "H:\\BetaGradient\\BetaScan\\BS529536\\Chiral\\BS02__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BSQuench\\Chiral\\BSQ__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BSQuench\\Chiral\\BSQ__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BGQuench\\Polyakov\\BGQ__{}_polyakov_ZSlice.csv"
# foldername = "I:\\BetaGradient\\Quenched\\BGQuench\\Polyakov\\BGQ__{}_polyakov_ZSlice.csv"
# foldername = "I:\\BetaGradient\\Quenched\\BSQuench\\Polyakov\\BSQ__{}_polyakov_ZSlice.csv"
foldername = "I:\\BetaGradient\\AtGradient\\AGQ03Polya\\AGQ__{}_polyakov_ZSlice.csv"
# lst = ["529", "530", "531", "532", "533", "534", "535", "536"]
lst = ["570", "575", "580", "585", "590", "595", "600", "605"]
savehead = "03"

polyalst = []
polyalste = []
susplst = []
susplste = []
times = []

polyalst2 = []
polyalste2 = []
susplst2 = []
susplste2 = []
times2 = []

for j in range(len(lst)):
    runname = "{}".format(lst[j])
    filename = foldername.format(lst[j])
    testarray = LoadMathematicaCSV(filename)
    arr = np.abs((testarray[:, 0] + testarray[:, 6]) / 2)
    v, s = JacknifeMean(arr)
    polyalst.append(v)
    polyalste.append(s)
    v, s = JacknifeCumulant(arr, "susp" + runname)
    susplst.append(v)
    susplste.append(s)
    v, s, t = AutoCorrelationSingleVariable(arr)
    times.append(t)
    arr = (np.abs(testarray[:, 0]) + np.abs(testarray[:, 6])) / 2
    v, s = JacknifeMean(arr)
    polyalst2.append(v)
    polyalste2.append(s)
    v, s = JacknifeCumulant(arr, "susp" + runname)
    susplst2.append(v)
    susplste2.append(s)
    v, s, t = AutoCorrelationSingleVariable(arr)
    times2.append(t)

polyalst = np.array(polyalst)
polyalste = np.array(polyalste)
susplst = np.array(susplst)
susplste = np.array(susplste)
times = np.array(times)

polyalst2 = np.array(polyalst2)
polyalste2 = np.array(polyalste2)
susplst2 = np.array(susplst2)
susplste2 = np.array(susplste2)
times2 = np.array(times2)


print(PrintAsMatlabArray(polyalst / 3, "polya" + savehead))
print(PrintAsMatlabArray(polyalste * np.sqrt(2 * times) / 3, "polyae" + savehead))
print(PrintAsMatlabArray(216 * 12 * 12 * 2 * susplst / 9, "polyasusp" + savehead))
print(PrintAsMatlabArray(216 * 12 * 12 * 2 * susplste * np.sqrt(2 * times) / 9, "polyasuspe" + savehead))
print(PrintAsMatlabArray(times, "ptime" + savehead))

print(PrintAsMatlabArray(polyalst2 / 3, "polya2" + savehead))
print(PrintAsMatlabArray(polyalste2 * np.sqrt(2 * times2) / 3, "polyae2" + savehead))
print(PrintAsMatlabArray(216 * 12 * 12 * 2 * susplst2 / 9, "polyasusp2" + savehead))
print(PrintAsMatlabArray(216 * 12 * 12 * 2 * susplste2 * np.sqrt(2 * times2) / 9, "polyasuspe2" + savehead))
print(PrintAsMatlabArray(times2, "ptime2" + savehead))
