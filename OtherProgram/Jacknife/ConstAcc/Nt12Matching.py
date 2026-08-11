import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeWilsonLoop, PrintAsMathematicaArray

folderHead = "H:\\ConstAccLarge\\Nt12\\Wilson{}\\"

testarrayr = LoadMathematicaCSV(folderHead.format(1) + "BSQ__VR_R.csv")
print(len(testarrayr))
print(testarrayr[46]*testarrayr[46])

# """
r0lst = []
r0elst = []
r1lst = []
r1elst = []
c0lst = []
c0elst = []
for i in range(9):
    for j in range(5):
        if 0 == j:
            testarray = LoadMathematicaCSV(folderHead.format(j + 1) + "BSQ__VR_Nt48_O{}.csv".format(i))
        else:
            testarray = np.vstack((testarray, LoadMathematicaCSV(folderHead.format(j + 1) + "BSQ__VR_Nt48_O{}.csv".format(i))))
    print(np.shape(testarray))
    r0, r0e, r1, r1e, c0, c0e = JacknifeWilsonLoop(testarray, testarrayr, 24, 2, 8, len(testarrayr), "O-{}".format(i), True)
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


