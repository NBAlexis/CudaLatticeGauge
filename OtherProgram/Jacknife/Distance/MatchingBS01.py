import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeWilsonLoop, PrintAsMathematicaArray

folderHead = "H:\\BS01\\Wilson\\"

testarrayr = LoadMathematicaCSV(folderHead + "BS01__VR_R.csv")
print(PrintAsMathematicaArray(testarrayr, "r"))
print(len(testarrayr))
print(testarrayr[44]*testarrayr[44])
betalst = ["53", "535", "54", "545", "55", "555", "56", "565", "57", "575", "58"]
tstart = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
tend =   [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
maxR =   [10, 10, 17, 17, 21, 21, 21, 21, 21, 21, 21]
# """
r0lst = []
r0elst = []
r1lst = []
r1elst = []
c0lst = []
c0elst = []
for i in range(0, 11):
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
r0={1.95074961,2.18662242,2.51266232,2.94907541,3.39833616,3.84148928,4.30994382,4.85433123,5.28256933,5.70878847,6.20695728};
r0e={0.01790624,0.01073408,0.00839091,0.00497807,0.00503037,0.0042662,0.00422978,0.00426579,0.00548036,0.00531822,0.00586454};
r1={1.2932804,1.50688747,1.76961799,2.08878002,2.41955867,2.74140142,3.07997736,3.47035468,3.78820828,4.10492205,4.47012746};
r1e={0.02468208,0.0038355,0.00524946,0.00170982,0.00178357,0.00204115,0.00233081,0.0025094,0.003503,0.00343051,0.00387983};
c0={-0.76148587,-0.71791597,-0.69922448,-0.71008132,-0.71002434,-0.70967176,-0.70901866,-0.70943487,-0.69891635,-0.68840307,-0.68096302};
c0e={0.04714158,0.01506971,0.01159282,0.00307568,0.00225402,0.0011375,0.00073049,0.00051554,0.00043288,0.00036241,0.00031202};
"""
