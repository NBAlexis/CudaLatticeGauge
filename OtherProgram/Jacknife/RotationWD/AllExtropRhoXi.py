"""
all used r-distributions

only for 502161 and 535161

only for polyakov, chiral, JG, JFL, JGS, JFS
only for omega = 0, 5, 10
"""
import numpy as np

from RotationWD.UsefulFunctions import MeasureValueType, getRDist, getRhoXi
from Visualization import errorbar, MarkerStyle, LineStyle

diskname = "G:\\RotationWD\\"
savefolder = "G:\\RotationWD\\ExportCSV\\RDistri\\"
nts = [3, 4, 5, 6]

fs25153 = ["Chiral"]
fs502161 = ["Chiral"]
fs5283161 = ["Chiral"]
fs535161 = ["Chiral1", "Chiral2"]
fsnames = [fs25153, fs502161, fs5283161, fs535161]

rhomeasuretypes = [
    MeasureValueType.JG,
    MeasureValueType.JFL
]

ximeasuretypes = [
    MeasureValueType.JGS,
    MeasureValueType.JFS
]

for i in range(4):
    vs = []
    ss = []
    for measuretype in rhomeasuretypes:
        v, s = getRhoXi(diskname, fsnames[i], measuretype, nts[i])
        vs.append(v)
        ss.append(s)
    vs = np.array(vs)
    ss = np.array(ss)
    np.savetxt(savefolder + "rho_nt{}.csv".format(nts[i]), np.vstack((vs, ss)), delimiter=',')

for i in range(4):
    vs = []
    ss = []
    for measuretype in ximeasuretypes:
        v, s = getRhoXi(diskname, fsnames[i], measuretype, nts[i])
        vs.append(v)
        ss.append(s)
    vs = np.array(vs)
    ss = np.array(ss)
    np.savetxt(savefolder + "xi_nt{}.csv".format(nts[i]), np.vstack((vs, ss)), delimiter=',')




