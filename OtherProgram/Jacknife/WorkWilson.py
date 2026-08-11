import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeWilsonLoop, PrintAsMathematicaArray

folderHead = "I:\\TreeImprove\\BSQ48\\Wilson48{}\\"
filelist = ["9", "10", "11"]
filelist2 = ["A", "B", "C", "D", "E"]

testarrayr = LoadMathematicaCSV(folderHead.format(filelist2[0]) + "BSQT__VR_R.csv")
print(len(testarrayr))
print(testarrayr[24]*testarrayr[24])
# """
r0lst = []
r0elst = []
r1lst = []
r1elst = []
c0lst = []
c0elst = []
for file in filelist:
    for i in range(5):
        if 0 == i:
            testarray = LoadMathematicaCSV(folderHead.format(filelist2[i]) + "BSQT__VR_Nt48_O{}.csv".format(file))
        else:
            testarray = np.vstack((testarray, LoadMathematicaCSV(folderHead.format(filelist2[i]) + "BSQT__VR_Nt48_O{}.csv".format(file))))
    print(np.shape(testarray))
    r0, r0e, r1, r1e, c0, c0e = JacknifeWilsonLoop(testarray, testarrayr, 24, 2, 8, len(testarrayr), file, True)
    print(PrintAsMathematicaArray([r0, r0e, r1, r1e, c0, c0e]))
    r0lst.append(r0)
    r0elst.append(r0e)
    r1lst.append(r1)
    r1elst.append(r1e)
    c0lst.append(c0)
    c0elst.append(c0e)

print(PrintAsMathematicaArray(r0lst))
print(PrintAsMathematicaArray(r0elst))
print(PrintAsMathematicaArray(r1lst))
print(PrintAsMathematicaArray(r1elst))
print(PrintAsMathematicaArray(c0lst))
print(PrintAsMathematicaArray(c0elst))
# """


