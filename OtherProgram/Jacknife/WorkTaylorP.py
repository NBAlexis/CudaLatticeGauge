from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads1 = "E:/TreeImprovedPolya/"
folderHeads2 = "E:/TreeImprovedPolya/"
# folderHeads = "E:/QuenchedResults/Nt48-56-65/GaugeMomentum/"
# folderHeads = "F:\\Builds\\FloatBuild-11-01\\Release\\"

lst = ["67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78"]
lst2 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
        "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "I2", "J2", "K2", "L2", "M2", "N2", "O2", "P2", "Q2", "R2", "S2", "T2", "U2", "V2", "W2", "X2", "Y2"]

values1 = []
stds1 = []
values2 = []
stds2 = []
values3 = []
stds3 = []
for i in range(len(lst)):
    for j in range(len(lst2)):
        fileNames1 = folderHeads1 + "Polyakov{}/BSQT__{}_polyakov.csv".format(lst2[j], lst[i])
        fileNames2 = folderHeads2 + "Angular{}/BSQT__angularJG_Nt6_All_O{}.csv".format(lst2[j], i)
        fileNames3 = folderHeads2 + "Angular{}/BSQT__angularJGS2_Nt6_All_O{}.csv".format(lst2[j], i)
        testarray1 = LoadMathematicaCSV(fileNames1)
        testarray2 = LoadMathematicaCSV(fileNames2)
        testarray3 = LoadMathematicaCSV(fileNames3)
        if 0 == j:
            testarrayre1 = np.abs(testarray1)
            testarrayre2 = np.real(testarray2) * 12 * 12 * 12 * 6
            testarrayre3 = np.real(testarray3) * 12 * 12 * 12 * 6
        else:
            testarrayre1 = np.hstack((testarrayre1, np.abs(testarray1)))
            testarrayre2 = np.hstack((testarrayre2, np.real(testarray2) * 12 * 12 * 12 * 6))
            testarrayre3 = np.hstack((testarrayre3, np.real(testarray3) * 12 * 12 * 12 * 6))
    print(len(testarrayre1))
    part1 = testarrayre1
    part2 = (-2*testarrayre3 + testarrayre2*testarrayre2)
    part3 = testarrayre1 * (-2*testarrayre3 + testarrayre2*testarrayre2)
    # v, s = JacknifeMean(part1, lst[i] + " P")
    v1 = np.mean(part1)
    s1 = np.std(part1) / np.sqrt(len(testarrayre1))
    print("==============", i)
    print(v1)
    print(s1)
    values1.append(v1)
    stds1.append(s1)
    print("==============", i)
    # v, s = JacknifeMean(part2, lst[i] + " S")
    v2 = np.mean(part2)
    s2 = np.std(part2) / np.sqrt(len(testarrayre1))
    print(v2)
    print(s2)
    values2.append(v2)
    stds2.append(s2)
    print("==============", i)
    # v, s = JacknifeMean(part3, lst[i] + " All")
    v3 = np.mean(part3)
    s3 = np.std(part3) / np.sqrt(len(testarrayre1))
    print(v3)
    print(s3)
    print(v1 * v2 - v3)
    values3.append(v3)
    stds3.append(s3)

print(PrintAsMathematicaArray(values1, "p"))
print(PrintAsMathematicaArray(stds1), "pe")
print(PrintAsMathematicaArray(values2, "d1"))
print(PrintAsMathematicaArray(stds2, "d1e"))
print(PrintAsMathematicaArray(values3, "d2"))
print(PrintAsMathematicaArray(stds3, "d2e"))

