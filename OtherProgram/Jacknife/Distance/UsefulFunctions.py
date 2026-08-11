import re
import numpy as np

def findK(folderStart, fileName):
    distlst = []
    with open(f'{folderStart}WilsonLoops\\{fileName}.txt', 'r') as f:
        for line_num, line in enumerate(f, start=1):
            if 0 != (line_num & 1):
                matches = re.findall(r"k\sis\s([\d]+)", line)
                distlst.append(int(matches[0]))
    return np.array(distlst)

def findDistances(folderStart, fileName):
    distlst = []
    with open(f'{folderStart}WilsonLoops\\{fileName}.txt', 'r') as f:
        for line_num, line in enumerate(f, start=1):
            if 0 != (line_num & 1):
                matches = re.findall(r"distance\sis\s([\d]+)", line)
                distlst.append(int(matches[0]))
    return np.array(distlst)
