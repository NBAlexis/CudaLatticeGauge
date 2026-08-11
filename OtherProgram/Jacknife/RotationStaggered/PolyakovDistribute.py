from RotationStaggered.UsefulFunctions import getRDist, MeasureValueType
from Visualization import errorbar

disk = "I"
names = ["534081162", "538078156", "542075015", "546072144"]
folders = ["Polyakov"]
nameidx = 1
name = names[nameidx]

pr, v, s = getRDist(disk, name, folders, MeasureValueType.Polya)
print(pr)
print(v)
print(s)

errorbar(pr, v, s)