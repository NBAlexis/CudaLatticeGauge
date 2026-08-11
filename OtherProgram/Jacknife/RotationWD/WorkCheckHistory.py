import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import JacknifeMean, JacknifeCumulant, LoadMathematicaCSV, PrintAsMatlabArray
from Visualization import errorbar, LineStyle, MarkerStyle, history

filenames = "F:\\Builds\\FloatBuild-09-11\\Release\\{}\\NC{}__polyakov_Nt3_{}.csv"
folders = ["47816"]
indexstart = 101
onlycheckerror = True

for f in folders:
    # start
    fout = filenames.format(f, f, "Out")
    dout = LoadMathematicaCSV(fout)
    errora = abs(dout) < 1.0e-6
    errorb = abs(dout) > 2.0
    errorc = np.isnan(dout)
    errord = np.isinf(dout)
    errorall = errora | errorb | errorc | errord
    print(fout)
    print(dout[np.where(errorall)])
    error1 = np.transpose(np.where(errorall))
    for i in range(20):
        print("errors:", i, (error1[error1[:, 0] == i][:, 1] + indexstart).tolist())
    print(len(dout[0]))
    if onlycheckerror:
        continue
    print(PrintAsMatlabArray(np.abs(dout[:, len(dout[0]) - 1])))
    for i in range(20):
        history([np.abs(dout[i, :])])
