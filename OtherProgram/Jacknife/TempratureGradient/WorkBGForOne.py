import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMathematicaArray

# foldername = "G:/BGQuench/Polyakov"
foldername = "G:/BGQuench/Chiral"
# lst1 = ["529", "530", "531", "532", "533", "534", "535", "536"]
lst1 = ["570", "575", "580", "585", "590", "595", "600", "605"]
# lst2 = ["Chiral1", "Chiral2"]
# lst2 = ["Polyakov1", "Polyakov2"]

v1lst = []
s1lst = []
v2lst = []
s2lst = []
for filename in lst1:
    f = 0
    # csvfile = foldername + "{}/BG02__{}_condensateZSlicepCCLightChiralKS.csv".format(lst2[f], filename)
    # csvfile = foldername + "{}/BSQ__{}_polyakov_ZSlice.csv".format(lst2[f], filename)
    # csvfile = foldername + "/BGQ__{}_polyakov_ZSlice.csv".format(filename)
    csvfile = foldername + "/BGQ__{}_condensateZSlicepCCLightChiralKS.csv".format(filename)

    testarray = LoadMathematicaCSV(csvfile)
    # f = 1
    # csvfile = foldername + "{}/BG02__{}_condensateZSlicepCCLightChiralKS.csv".format(lst2[f], filename)
    # csvfile = foldername + "{}/BSQ__{}_polyakov_ZSlice.csv".format(lst2[f], filename)
    # testarray = np.vstack((testarray, LoadMathematicaCSV(csvfile) * -1))
    # testarray = np.vstack((testarray, LoadMathematicaCSV(csvfile)))
    position0 = -np.real(testarray[:, 0])
    position6 = -np.real(testarray[:, 6])
    v1, s1 = JacknifeMean(position0)
    v2, s2 = JacknifeMean(position6, filename)
    v1lst.append(v1)
    v2lst.append(v2)
    s1lst.append(s1)
    s2lst.append(s2)

print(PrintAsMathematicaArray(v1lst, "polya021"))
print(PrintAsMathematicaArray(s1lst, "polya021e"))
print(PrintAsMathematicaArray(v2lst, "polya027"))
print(PrintAsMathematicaArray(s2lst, "polya027e"))




