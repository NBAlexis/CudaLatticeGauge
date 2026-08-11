import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMatlabArray
from RotationWD.UsefulFunctions import getMean, getPrAndNumberOfSites, MeasureValueType
from Visualization import errorbar


diskname = "G:\\RotationWD\\"
savefolder = "G:\\RotationWD\\ExportCSV\\Extrop\\"
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

# signs are reductent, we have used new measurement
signs25153 = [1]
signs502161 = [1]
signs5283161 = [1]
signs535161 = [1, 1]
signs = [signs25153, signs502161, signs5283161, signs535161]

for i in range(4):
    for measuretype in MeasureValueType:
        if measuretype == MeasureValueType.Polya:
            print(measuretype, i)
            getMean(diskname, pfsnames[i], [1], measuretype, nts[i], savefolder)
        else:
            print(measuretype, i)
            getMean(diskname, fsnames[i], signs[i], measuretype, nts[i], savefolder)

