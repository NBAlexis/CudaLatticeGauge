from enum import Enum

import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean


class MeasureValueType(Enum):
    Polya = 0
    ChiralUD = 1
    ChiralS = 2
    JG = 3
    JGChen = 4
    JGS = 5
    JGPot = 6
    JGSurf = 7
    JFLUD = 8
    JFLS = 9
    JFSUD = 10
    JFSS = 11
    JFPotUD = 12
    JFPotS = 13

def getFileName(measureVal: MeasureValueType) -> str:
    if measureVal is MeasureValueType.Polya:
        return "polyakov_OverR"
    elif measureVal is MeasureValueType.ChiralUD:
        return "condensatepCCHeavyChiralKS_OverR"
    elif measureVal is MeasureValueType.ChiralS:
        return "condensatepCCLightChiralKS_OverR"
    elif measureVal is MeasureValueType.JG:
        return "angularJG"
    elif measureVal is MeasureValueType.JGChen:
        return "angularJGChen"
    elif measureVal is MeasureValueType.JGS:
        return "angularJGS"
    elif measureVal is MeasureValueType.JGPot:
        return "angularJGPot"
    elif measureVal is MeasureValueType.JGSurf:
        return "angularJGSurf"
    elif measureVal is MeasureValueType.JFLUD:
        return "condensatepFAHeavyOrbitalKS"
    elif measureVal is MeasureValueType.JFLS:
        return "condensatepFALightOrbitalKS"
    elif measureVal is MeasureValueType.JFSUD:
        return "condensatepFAHeavySpinKS"
    elif measureVal is MeasureValueType.JFSS:
        return "condensatepFALightSpinKS"
    elif measureVal is MeasureValueType.JFPotUD:
        return "condensatepFAHeavyPotentialKS"
    elif measureVal is MeasureValueType.JFPotS:
        return "condensatepFALightPotentialKS"

"""
def getSaveName(measureVal: MeasureValueType) -> str:
    if measureVal is MeasureValueType.Polya:
        return "polya"
    elif measureVal is MeasureValueType.Chiral:
        return "chiral"
    elif measureVal is MeasureValueType.JG:
        return "JG"
    elif measureVal is MeasureValueType.JGChen:
        return "JGChen"
    elif measureVal is MeasureValueType.JGS:
        return "JGS"
    elif measureVal is MeasureValueType.JGPot:
        return "JGPot"
    elif measureVal is MeasureValueType.JGSurf:
        return "JGSurf"
    elif measureVal is MeasureValueType.JFL:
        return "JFL"
    elif measureVal is MeasureValueType.JFS:
        return "JFS"
    elif measureVal is MeasureValueType.JFPot:
        return "JFPot"
"""

def getPrAndNumberOfSites(Lx: int):
    """

    :param Lx:
    :return:
    pr2: real radious for each with R data
    numSite2:
    np.sum(numSite2): volume
    idxedge:
    """
    nEdge = (Lx - 1) # 11
    halfX = Lx // 2
    numSite1 = np.array([0 for _ in range(4 * Lx * Lx + 1)])
    pr1 = np.array([0 for _ in range(4 * Lx * Lx + 1)])
    for x in range(0, Lx):
        for y in range(0, Lx):
            xx = 2 * (x - halfX) + 1
            yy = 2 * (y - halfX) + 1
            prSq = xx * xx + yy * yy
            pr1[prSq] = prSq
            numSite1[prSq] = numSite1[prSq] + 1
    pr2 = pr1[numSite1 > 0]
    numSite2 = numSite1[numSite1 > 0]
    idxedge = np.where(pr2 == nEdge * nEdge + 1)[0][0] - 1
    pr2 = np.sqrt(pr2) / 2
    return pr2, numSite2, np.sum(numSite2), idxedge

def getFileNamesByType(diskName, header, folderName, omega: int, measureVal: MeasureValueType):
    if measureVal is MeasureValueType.Polya or measureVal is MeasureValueType.ChiralUD or measureVal is MeasureValueType.ChiralS:
        return f"{diskName}:\\NewSRF\\SRF{header}\\{folderName}\\SRF{header}__{omega}_{getFileName(measureVal)}.csv"
    else:
        return f"{diskName}:\\NewSRF\\SRF{header}\\{folderName}\\SRF{header}__{getFileName(measureVal)}_Nt6_O{omega}.csv"

def getRDist(diskName, header, folderNames, measureVal: MeasureValueType):
    """
    Polyakov use abs
    """
    polyakovUseAbsOrReal = True
    resvall = []
    ressall = []
    pr, nums, _, _ = getPrAndNumberOfSites(12)
    for omega in range(0, 16):
        # if measureVal is not MeasureValueType.Polya and measureVal is not MeasureValueType.Chiral:
        #     if 0 == omega:
        #         continue
        resv = []
        ress = []
        arr = LoadMathematicaCSV(getFileNamesByType(diskName, header, folderNames[0], omega, measureVal))
        if len(folderNames) > 1:
            for j in range(1, len(folderNames)):
                arrj = LoadMathematicaCSV(getFileNamesByType(diskName, header, folderNames[j], omega, measureVal))
                arr = np.vstack((arr, arrj))
        if measureVal is MeasureValueType.Polya:
            rearAll = np.abs(arr) if polyakovUseAbsOrReal else np.real(arr)
            rearAll = rearAll / 3
        elif measureVal is MeasureValueType.ChiralUD or measureVal is MeasureValueType.ChiralS:
            rearAll = np.real(arr)
        else:
            rearAll = np.real(arr)
        for i in range(len(pr)):
            v, s = JacknifeMean(rearAll[:, i])
            _, _, t = AutoCorrelationSingleVariable(rearAll[:, i])
            resv.append(v)
            ress.append(s * np.sqrt(2 * t))
        resvall.append(resv)
        ressall.append(ress)
    return pr, resvall, ressall