import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import JacknifeMean, JacknifeCumulant, LoadMathematicaCSV, PrintAsMatlabArray
from Visualization import errorbar, LineStyle, MarkerStyle

filenames = "H:\\configs\\newresf\\{}\\Polyakov\\NC{}__polyakov_Nt3_{}.csv"
folders = ["49916"]

for f in folders:
    # start
    fout = filenames.format(f, f, "Out")
    dout = LoadMathematicaCSV(fout)
    print(len(dout[0]))
    print(PrintAsMatlabArray(np.abs(dout[:, len(dout[0]) - 1])))

