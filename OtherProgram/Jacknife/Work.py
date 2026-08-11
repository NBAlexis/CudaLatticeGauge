from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

# folderHeads = "E:/QuenchedResults/Nt6-56-65/Polyakov/"
# folderHeads = "E:/QuenchedResults/Nt6-56-65/GaugeMomentum/"
# folderHeads = "E:/QuenchedResults/Nt48-56-65/GaugeMomentum/"
# folderHeads = "H:\\BS2\\529\\"
# folderHeads = "H:\\BS2\\536\\"
# folderHeads = "F:\\Builds\\FloatBuild-11-01\\Release\\"
# folderHeads = "G:\\BetaGradient\\Quenched\\BSQuench\\Polyakov\\"
folderHeads = "G:\\MesonMassElectric\\BS01\\Polyakov\\"
# folderHeads = "G:\\MesonMassElectric\\BS01\\Chiral\\"
# folderHeads = "E:\\QuenchedResults\\BSQT\\"

# lst = ["56", "57", "58", "59", "60", "61", "62", "63", "64", "65"]
lst = ["52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63"]
# lst = ["67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78"]
# lst = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
# lst = ["570", "575", "580", "585", "590", "595", "600", "605"]
# lst2 = ["BSQ", "BSQ2", "BSQ3", "BSQ4", "BSQ5", "BSQ6", "BSQ7"]
# lst = ["51400", "51775", "51900", "52150", "52400", "52525", "52900", "53275", "53400", "53650", "53900", "54025", "54400"]
# lst = ["52100", "52475", "52600", "52850", "53100", "53225", "53600", "53975", "54100", "54350", "54600", "54725", "55100"]
# lst2 = ["0-30", "30-60"]
lst2 = [""]
# lst2 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "T", "S", "T", "U"]

values = []
stds = []
tvalues = []
for i in range(len(lst)):
    for j in range(len(lst2)):
        # fileNames = folderHeads + "{}__{}_polyakov.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJG_Nt6_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJGS2_Nt6_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJG_Nt48_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJGS2_Nt48_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "Polyakov\\BS529__{}_polyakov.csv".format(lst[i], lst[i])
        # fileNames = folderHeads + "Polyakov\\BS536__{}_polyakov.csv".format(lst[i], lst[i])
        fileNames = folderHeads + "BS01__{}_polyakov.csv".format(lst[i])
        # fileNames = folderHeads + "BS01__{}_condensatepCCLightChiralKS.csv".format(lst[i])
        # fileNames = folderHeads + "Polyakov{}\\BSQT__{}_polyakov.csv".format(lst2[j], lst[i])

        testarray = LoadMathematicaCSV(fileNames)
        if 0 == j:
            testarrayre = np.abs(testarray)
            # testarrayre = np.real(testarray)
        else:
            testarrayre = np.hstack((testarrayre, np.abs(testarray)))
            # testarrayre = np.hstack((testarrayre, np.real(testarray)))
    # testarrayre = testarrayre * testarrayre
    # v, s = JacknifeMean(testarrayre, lst[i])
    v, s = JacknifeCumulant(testarrayre, lst[i])
    # _, _, t = AutoCorrelationSingleVariable(testarrayre)
    print("==============", i, "with {} results".format(len(testarrayre)))
    print(v)
    print(s)
    # print(t)
    values.append(v)
    stds.append(s)
    # tvalues.append(t)


print(PrintAsMathematicaArray(values, "susp"))
print(PrintAsMathematicaArray(stds, "suspe"))
# print(PrintAsMathematicaArray(tvalues, "t"))

