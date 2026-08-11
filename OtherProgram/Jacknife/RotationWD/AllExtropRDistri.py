"""
all used r-distributions

only for 502161 and 535161

only for polyakov, chiral, JG, JFL, JGS, JFS
only for omega = 0, 5, 10
"""
from RotationWD.UsefulFunctions import MeasureValueType, getRDist
from Visualization import errorbar, MarkerStyle, LineStyle

diskname = "G:\\RotationWD\\"
savefolder = "G:\\RotationWD\\ExportCSV\\RDistri\\"
nts = [3, 4, 5, 6]

fs25153 = ["Chiral"]
fs502161 = ["Chiral"]
fs5283161 = ["Chiral"]
fs535161 = ["Chiral1", "Chiral2"]
fsnames = [fs25153, fs502161, fs5283161, fs535161]

pfs25153 = ["Polyakov"]
pfs502161 = ["Polyakov"]
pfs5283161 = ["Polyakov"]
pfs535161 = ["Polyakov"]
pfsnames = [pfs25153, pfs502161, pfs5283161, pfs535161]

measuretypes = [
    MeasureValueType.Polya,
    MeasureValueType.Chiral,
    MeasureValueType.JG,
    MeasureValueType.JFL,
    MeasureValueType.JGS,
    MeasureValueType.JFS
]

omegalst = [0, 5, 10]

for i in range(4):
    for measuretype in measuretypes:
        if measuretype == MeasureValueType.Polya:
            getRDist(diskname, pfsnames[i], measuretype, nts[i], omegalst, savefolder)
        else:
            getRDist(diskname, fsnames[i], measuretype, nts[i], omegalst, savefolder)


