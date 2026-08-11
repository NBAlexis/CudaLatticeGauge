import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeWilsonLoop, PrintAsMathematicaArray

folderHead = "H:\\Quench242424-48\\Wilson\\"

testarrayr = LoadMathematicaCSV(folderHead + "BSQ__VR_R.csv")
print(PrintAsMathematicaArray(testarrayr, "r"))
print(len(testarrayr))
print(testarrayr[44]*testarrayr[44])
betalst = ["565", "57", "575", "58", "585", "59", "595", "60", "605", "61", "615"]
tstart = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
tend =   [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
maxR =   [21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21]
# """
r0lst = []
r0elst = []
r1lst = []
r1elst = []
c0lst = []
c0elst = []
for i in range(0, 11):
    testarray = LoadMathematicaCSV(folderHead + "BSQ__VR_Nt48_{}.csv".format(betalst[i]))
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
r0={2.58745983,2.93412439,3.30216761,3.68752283,4.09156952,4.49300951,4.93543565,5.32649375,5.81808429,6.24351864,6.64796299};
r0e={0.00948706,0.00625953,0.00472433,0.0038928,0.0036356,0.00396548,0.00417737,0.00431635,0.00480782,0.00561304,0.00569953};
r1={1.86758342,2.10182481,2.36923367,2.64656421,2.9426458,3.23900849,3.56002773,3.84948111,4.21078519,4.52621793,4.82836829};
r1e={0.01054109,0.00211354,0.00158517,0.00176099,0.00191489,0.00235542,0.00257202,0.00277809,0.00315774,0.00376024,0.00385669};
c0={-0.62865279,-0.66071031,-0.66709412,-0.67117189,-0.66914038,-0.66394305,-0.66186516,-0.65403023,-0.6483903,-0.640316,-0.6311145};
c0e={0.01987097,0.00557506,0.00241499,0.00123629,0.00079085,0.0005626,0.00044967,0.00035942,0.00031121,0.00028277,0.0002421};
"""
