"""
we do only for Polyakov, susp of Polyakov, chiral and susp of chiral
"""
import numpy as np

from RotationWD.UsefulFunctions import MeasureValueType, getVDist
from Visualization import errorbar, LineStyle, MarkerStyle

diskname = "G:\\RotationWD\\"
savefolder = "G:\\RotationWD\\ExportCSV\\Extrop\\"
nts = [3, 4, 5, 6]
edges = [5, 7, 9, 11]

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

for i in range(4):
  x, v, s, sv, ss = getVDist(diskname, pfsnames[i], MeasureValueType.Polya, nts[i], edges[i])
  _, v2, s2, sv2, ss2 = getVDist(diskname, fsnames[i], MeasureValueType.Chiral, nts[i], edges[i])
  filetosave = savefolder + "velocity_nt" + str(nts[i]) + ".csv"
  print("saving:", filetosave)
  np.savetxt(filetosave, np.vstack((
      np.array(x),
      np.array(v),
      np.array(s),
      np.array(sv),
      np.array(ss),
      np.array(v2),
      np.array(s2),
      np.array(sv2),
      np.array(ss2))), delimiter=',')


