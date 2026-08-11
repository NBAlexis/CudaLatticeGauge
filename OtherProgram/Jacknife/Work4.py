import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, PrintAsMathematicaArray

folderHeads = "E:/QuenchedResults/QuenchRotation/56QT/Polyakov/"

lst1 = ["A", "B"]
lst2 = ["", "2"]

testarrayCombine = None
for i in range(len(lst1)):
    fileNames = folderHeads + "{}/56QT{}__polyakov_Nt6_In.csv".format(lst1[i], lst2[i])
    testarray = LoadMathematicaCSV(fileNames)
    if 0 == i:
        testarrayCombine = testarray
    else:
        testarrayCombine = np.hstack((testarrayCombine, testarray))

values = []
stds = []
for i in range(len(testarrayCombine)):
    v, s = JacknifeMean(np.abs(testarrayCombine[i,:]))
    print("==============", i)
    print(v)
    print(s)
    values.append(v)
    stds.append(s)

print(PrintAsMathematicaArray(np.array(values)/3))
print(PrintAsMathematicaArray(np.array(stds)/3))
