from io import StringIO

import numpy as np


def LoadMathematicaCSV(fileName: str):
    with open(fileName) as f:
        contents = f.read()
        contents = contents.replace("*^", "e")
        contents = contents.replace(" I", "j")
        contents = contents.replace(" ", "")
        c = StringIO(contents)
        d = np.loadtxt(c, delimiter=",", dtype=str)
        d = d.astype(complex)
        return d


def JacknifeMean(lst):
    size = len(lst)
    orignalv = np.mean(lst)
    argsv = 0
    avers = []
    for i in range(size):
        subseto = np.delete(lst, i)
        avers.append(np.mean(subseto))
        argsv = argsv + (avers[i] - orignalv) ** 2
    argsv = argsv * (size - 1) / size
    unbaisv = orignalv + (size - 1) * (orignalv - np.mean(np.array(avers)))
    return unbaisv, np.sqrt(argsv)

def JacknifeCumulant(lst):
    size = len(lst)
    orignalv = np.mean(lst * lst) - np.mean(lst) * np.mean(lst)
    argsv = 0
    avers = []
    for i in range(size):
        subseto = np.delete(lst, i)
        avers.append(np.mean(subseto * subseto) - np.mean(subseto) * np.mean(subseto))
        argsv = argsv + (avers[i] - orignalv) ** 2
    argsv = argsv * (size - 1) / size
    unbaisv = orignalv + (size - 1) * (orignalv - np.mean(np.array(avers)))
    return unbaisv, np.sqrt(argsv)
