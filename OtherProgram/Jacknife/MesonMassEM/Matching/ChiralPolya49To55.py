from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

folderHeads = "H:\\BS01Chiral1\\"
# folderHeads = "H:\\BS01Polyakov\\"


lst = ["49", "50", "51", "52", "53", "54", "55", "56", "57"]
# lst = ["67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78"]
# lst = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
# lst = ["570", "575", "580", "585", "590", "595", "600", "605"]
# lst2 = ["BSQ", "BSQ2", "BSQ3", "BSQ4", "BSQ5", "BSQ6", "BSQ7"]
# lst = ["51400", "51775", "51900", "52150", "52400", "52525", "52900", "53275", "53400", "53650", "53900", "54025", "54400"]
# lst = ["52100", "52475", "52600", "52850", "53100", "53225", "53600", "53975", "54100", "54350", "54600", "54725", "55100"]
# lst2 = ["0-30", "30-60"]
lst2 = [""]
# lst2 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "T", "S", "T", "U"]

cv = []
sv = []
values = []
stds = []
tvalues = []
for i in range(len(lst)):
    for j in range(len(lst2)):
        # fileNames = folderHeads + "BS01__{}_polyakov.csv".format(lst[i])
        fileNames = folderHeads + "BS01__{}_condensatepCCLightChiralKS.csv".format(lst[i])

        testarray = LoadMathematicaCSV(fileNames)
        if 0 == j:
            # testarrayre = np.abs(testarray)
            testarrayre = np.real(testarray)
        else:
            # testarrayre = np.hstack((testarrayre, np.abs(testarray)))
            testarrayre = np.hstack((testarrayre, np.real(testarray)))
    # testarrayre = testarrayre * testarrayre
    v, s = JacknifeMean(testarrayre, lst[i])
    _, _, t = AutoCorrelationSingleVariable(testarrayre)
    cv.append(-v)
    sv.append(s * np.sqrt(2 * t))
    v, s = JacknifeCumulant(testarrayre, lst[i])
    print("==============", i, "with {} results".format(len(testarrayre)))
    print(v)
    print(s)
    # print(t)
    values.append(v)
    stds.append(s * np.sqrt(2 * t))
    # tvalues.append(t)


print(PrintAsMathematicaArray(values, "susp"))
print(PrintAsMathematicaArray(stds, "suspe"))
# print(PrintAsMathematicaArray(tvalues, "t"))

errorbar(range(len(values)), [cv], [sv])
errorbar(range(len(values)), [values], [stds])

