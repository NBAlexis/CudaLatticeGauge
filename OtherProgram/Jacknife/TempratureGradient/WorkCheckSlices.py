import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMathematicaArray

# foldername = "H:\\BGQuench\\Polyakov\\BGQ__{}_polyakov_ZSlice.csv"
# foldername = "H:\\BGQuench\\Chiral\\BGQ__{}_condensateZSlicepCCLightChiralKS.csv"
# lst = ["570", "575", "580", "585", "590", "595", "600", "605"]

# foldername = "H:\\BetaGradient\\BetaScan\\BS529536\\Polyakov\\BS02__{}_polyakov_ZSlice.csv"
# foldername = "H:\\BetaGradient\\BetaScan\\BS529536\\Chiral\\BS02__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BSQuench\\Chiral\\BSQ__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BSQuench\\Chiral\\BSQ__{}_condensateZSlicepCCLightChiralKS.csv"
# foldername = "G:\\BGQuench\\Polyakov\\BGQ__{}_polyakov_ZSlice.csv"
foldername = "I:\\BetaGradient\\AtGradient\\AGQ02Polya\\AGQ__{}_polyakov_ZSlice.csv"
# lst = ["529", "530", "531", "532", "533", "534", "535", "536"]
lst = ["570", "575", "580", "585", "590", "595", "600", "605"]

polyalst = []

for j in range(len(lst)):
    runname = "{}".format(lst[j])
    filename = foldername.format(lst[j])
    testarray = LoadMathematicaCSV(filename)
    thisbetap = []
    for k in range(12):
        arr = np.abs(testarray[:, k])
        v, _ = JacknifeMean(arr, "{}-{}".format(lst[j], k))
        thisbetap.append(v)
    polyalst.append(thisbetap)
# print(PrintAsMathematicaArray(polyalst, "chiral"))
# print(PrintAsMathematicaArray(polyalste, "chirale"))
# print(PrintAsMathematicaArray(susplst, "chiralsusp"))
# print(PrintAsMathematicaArray(susplste, "chiralsuspe"))
# print(PrintAsMathematicaArray(times, "time"))
print(PrintAsMathematicaArray(polyalst, "polya"))
