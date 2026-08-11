from JacknifePrograms import *

folderHeads = "E:/QuenchedResults/Nt6-56-65/Polyakov/"
# folderHeads = "E:/QuenchedResults/Nt6-56-65/GaugeMomentum/"
# folderHeads = "E:/QuenchedResults/Nt48-56-65/GaugeMomentum/"
# folderHeads = "H:\\BS2\\529\\"
# folderHeads = "H:\\BS2\\536\\"

lst = ["56", "57", "58", "59", "60", "61", "62", "63", "64", "65"]
# lst = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
lst2 = ["BSQ", "BSQ2", "BSQ3", "BSQ4", "BSQ5", "BSQ6"]
# lst = ["51400", "51775", "51900", "52150", "52400", "52525", "52900", "53275", "53400", "53650", "53900", "54025", "54400"]
# lst = ["52100", "52475", "52600", "52850", "53100", "53225", "53600", "53975", "54100", "54350", "54600", "54725", "55100"]

values = []
stds = []
for i in range(len(lst)):
    for j in range(len(lst2)):
        fileNames = folderHeads + "{}__{}_polyakov.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJG_Nt6_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJGS2_Nt6_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJG_Nt48_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "{}__angularJGS2_Nt48_All_O{}.csv".format(lst2[j], lst[i])
        # fileNames = folderHeads + "Polyakov\\BS529__{}_polyakov.csv".format(lst[i], lst[i])
        # fileNames = folderHeads + "Polyakov\\BS536__{}_polyakov.csv".format(lst[i], lst[i])

        testarray = LoadMathematicaCSV(fileNames)
        if 0 == j:
            # testarrayre = np.real(testarray)
            testarrayre = np.abs(testarray)
        else:
            # testarrayre = np.hstack((testarrayre, np.real(testarray)))
            testarrayre = np.hstack((testarrayre, np.abs(testarray)))
    # testarrayre = testarrayre * testarrayre
    v, s = JacknifeMean(testarrayre, lst[i])
    # v, s = JacknifeCumulant(testarrayre)
    print("==============", i)
    print(v)
    print(s)
    values.append(v)
    stds.append(s)


print(PrintAsMathematicaArray(values))
print(PrintAsMathematicaArray(stds))
