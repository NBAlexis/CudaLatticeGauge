import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeWilsonLoop, PrintAsMathematicaArray

folderHead = "H:\\BS01Wilson\\"

testarrayr = LoadMathematicaCSV(folderHead + "BS01__VR_R.csv")
print(PrintAsMathematicaArray(testarrayr, "r"))
print(len(testarrayr))
print(testarrayr[49]*testarrayr[49])
betalst = ["49", "50", "51", "52", "53", "54", "55", "56", "57"]
tstart = [1, 1, 1, 2, 2, 2, 2, 2, 2]
tend =   [3, 3, 3, 6, 6, 8, 8, 8, 8]
maxR =   [6, 6, 6, 6, 8, 10, 15, 15, 15]
# """
r0lst = []
r0elst = []
r1lst = []
r1elst = []
c0lst = []
c0elst = []
for i in range(5, 6):
    testarray = LoadMathematicaCSV(folderHead + "BS01__VR_Nt48_{}.csv".format(betalst[i]))
    # testarray = np.abs(testarray[100:,:])
    # testarray = testarray[100:, :]
    print(np.shape(testarray))
    # print(testarray)
    r0, r0e, r1, r1e, c0, c0e = JacknifeWilsonLoop(testarray, testarrayr, 24, tstart[i], tend[i], maxR[i], betalst[i], True)
    print(PrintAsMathematicaArray([r0, r0e, r1, r1e, c0, c0e]))
    r0lst.append(r0)
    r0elst.append(r0e)
    r1lst.append(r1)
    r1elst.append(r1e)
    c0lst.append(c0)
    c0elst.append(c0e)
    print(r0 / r1)

print(PrintAsMathematicaArray(r0lst, "r0"))
print(PrintAsMathematicaArray(r0elst, "r0e"))
print(PrintAsMathematicaArray(r1lst, "r1"))
print(PrintAsMathematicaArray(r1elst, "r1e"))
print(PrintAsMathematicaArray(c0lst, "c0"))
print(PrintAsMathematicaArray(c0elst, "c0e"))
# """

"""

"""
