import numpy as np

from JacknifePrograms import LoadMathematicaCSV, JacknifeMean

# foldername = "I:/BetaGradient/BetaScan/BS529536/Polyakov/"
foldername = "I:/BetaGradient/BetaScan/BS529536/Polyakov/"
# lst = ["529", "530", "531", "532", "533", "534", "535", "536"]
lst = ["529", "530", "531", "532", "533", "534", "535", "536"]

for i in range(len(lst)):
    fileNames = foldername + "BS02__{}_polyakov.csv".format(lst[i])
    testarray = LoadMathematicaCSV(fileNames)
    v, s = JacknifeMean(np.real(testarray))
    print("==============", i)
    print(v)
    print(s)
    v, s = JacknifeMean(np.abs(testarray))
    print(v)
    print(s)
