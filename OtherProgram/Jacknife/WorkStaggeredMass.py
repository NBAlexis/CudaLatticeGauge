from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import *

disk = "H"
beta = ["490", "550"]
tobechecked1 = ["PS", "PV", "PV", "PV", "S", "VT", "VT", "VT"]
tobechecked2 = ["uu", "dd", "ud", "du"]
tobechecked3 = ["", "xy", "xz", "yz", "", "x", "y", "z"]
T = 12
d = 2

for dd in range(3):
    d = dd + 1
    for b in beta:
        allem = []
        for em in range(12):
            allud = []
            for ud in range(4):
                allpv = []
                for pv in range(8):
                    filename = "{}:\\MesonMassElectric\\Meson\\{}01\\EC{}\\EC{}_01__{}{}{}_{}.csv".\
                        format(disk, b, em, b, tobechecked1[pv], tobechecked2[ud], tobechecked3[pv], em)
                    testarray = LoadMathematicaCSV(filename)
                    testarray = testarray[100:1000, :]
                    res = JacknifeMeasonSimple(testarray, T, d)
                    print("{} {} {}{} {} = {}".format(b, em, tobechecked1[pv], tobechecked3[pv], tobechecked2[ud], res))
                    allpv.append(res)
                allud.append(allpv)
            allem.append(allud)
        alldata = np.reshape(np.array(allem), (12, 4 * 8 * 2))
        np.savetxt("{}:\\MesonMassElectric\\Meson\\mass{}-d{}.txt".format(disk, b, d), alldata, delimiter=',')

