"""
all condensation is measured as tr[D^-1], not multiplied kappa, not multiplied the minus sign

For gauge angular momentum, the beta/Nc is already multiplied
For fermion angular momentum, the -2kappa is multiplied, minus sign comes from -tr[D], 2 comes from Nf, and kappa is in definition.
"""
from enum import Enum

import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant


class MeasureValueType(Enum):
    Polya = 0
    Chiral = 1
    JG = 2
    JGChen = 3
    JGS = 4
    JGPot = 5
    JGSurf = 6
    JFL = 7
    JFS = 8
    JFPot = 9

def getPrAndNumberOfSites(Lx: int):
    nEdge = (Lx - 2) // 2
    numSite1 = np.array([0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    pr1 = np.array([0.0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    for x in range(-nEdge, nEdge + 1):
        for y in range(-nEdge, nEdge + 1):
            prSq = x * x + y * y
            pr1[prSq] = prSq
            numSite1[prSq] = numSite1[prSq] + 1
    pr2 = pr1[numSite1 > 0]
    numSite2 = numSite1[numSite1 > 0]
    idxedge = np.where(pr2 == nEdge * nEdge)[0][0]
    pr2 = np.sqrt(pr2)
    return pr2, numSite2, np.sum(numSite2), idxedge

def getPrAndNumberOfSitesTighter(Lx: int, edge: int):
    nEdge = (Lx - 2) // 2
    numSite1 = np.array([0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    pr1 = np.array([0.0 for _ in range(nEdge * nEdge + nEdge * nEdge + 1)])
    for x in range(-nEdge, nEdge + 1):
        for y in range(-nEdge, nEdge + 1):
            prSq = x * x + y * y
            pr1[prSq] = prSq
            numSite1[prSq] = numSite1[prSq] + 1
    pr2 = pr1[numSite1 > 0]
    numSite2 = numSite1[numSite1 > 0]
    idxedge = np.where(pr2 == edge * edge)[0][0]
    pr2 = np.sqrt(pr2)
    return pr2, numSite2, np.sum(numSite2), idxedge

def getFileName(measureVal: MeasureValueType) -> str:
    if measureVal is MeasureValueType.Polya:
        return "polyakov"
    elif measureVal is MeasureValueType.Chiral:
        return "condensateChiral"
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
    elif measureVal is MeasureValueType.JFL:
        return "angularJL"
    elif measureVal is MeasureValueType.JFS:
        return "angularJS"
    elif measureVal is MeasureValueType.JFPot:
        return "angularJPot"

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

def getHeadByNt(nt: int)-> str:
    if 3 == nt:
        return "25153"
    elif 4 == nt:
        return "502161"
    elif 5 == nt:
        return "5283161"
    elif 6 == nt:
        return "535161"

def getInVolume(nt: int) -> int:
    """
    Lx * Ly * Lz remove Dirichlet sites
    """
    if 3 == nt:
        return 654
    elif 4 == nt:
        return 1544
    elif 5 == nt:
        return 3050
    elif 6 == nt:
        return 5244

def getOutVolume(nt: int) -> int:
    """
    Lx * Ly * Lz
    """
    return 2 * nt * (4 * nt + 1) * (4 * nt + 1)

def getKappa(nt: int) -> float:
    if 3 == nt:
        return 0.153
    elif 4 == nt:
        return 0.161
    elif 5 == nt:
        return 0.161
    elif 6 == nt:
        return 0.161

def getBeta(nt: int) -> float:
    if 3 == nt:
        return 2.5
    elif 4 == nt:
        return 5.02
    elif 5 == nt:
        return 5.283
    elif 6 == nt:
        return 5.35

def getOmegaSep(nt: int) -> float:
    if 3 == nt:
        return 0.0117
    elif 4 == nt:
        return 0.008775
    elif 5 == nt:
        return 0.00702
    elif 6 == nt:
        return 0.00585

def getMean(diskName, folderNames, signs, measureVal: MeasureValueType, nt: int, savefolder: str):
    """
    For Polyakov, in, and out is used, abs is used. Real is used
    For others, radial distribution is used, and in/out is re-calculated.
    For chiral 2kappa is multiplied.
    For fermion angular momenta, Nf is divided.
    """
    polyakovUseAbsOrReal = False
    if measureVal is MeasureValueType.Polya:
        arrIn = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_In.csv".format(
              diskName,
              getHeadByNt(nt),
              nt,
              folderNames[0],
              getHeadByNt(nt),
              getFileName(measureVal),
              nt))
        arrOut = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_Out.csv".format(
               diskName,
               getHeadByNt(nt),
               nt,
               folderNames[0],
               getHeadByNt(nt),
               getFileName(measureVal),
               nt))
        arrIn = np.abs(arrIn) if polyakovUseAbsOrReal else np.real(arrIn)
        arrOut = np.abs(arrOut) if polyakovUseAbsOrReal else np.real(arrOut)
        polyaIn = []
        polyaIne = []
        suspIn = []
        suspIne = []
        polyaOut = []
        polyaOute = []
        suspOut = []
        suspOute = []
        for i in range(0, 11):
            v, s = JacknifeMean(arrIn[i])
            _, _, t = AutoCorrelationSingleVariable(arrIn[i])
            polyaIn.append(v)
            polyaIne.append(s * np.sqrt(2 * t))
            v, s = JacknifeCumulant(arrIn[i])
            suspIn.append(v)
            suspIne.append(s * np.sqrt(2 * t))
            v, s = JacknifeMean(arrOut[i])
            _, _, t = AutoCorrelationSingleVariable(arrOut[i])
            polyaOut.append(v)
            polyaOute.append(s * np.sqrt(2 * t))
            sv, ss = JacknifeCumulant(arrOut[i])
            suspOut.append(sv)
            suspOute.append(ss * np.sqrt(2 * t))
        outv = getOutVolume(nt) / 9
        inv = getInVolume(nt) / 9
        savefilename = savefolder + getSaveName(measureVal) + "_nt{}.csv".format(nt)
        print(savefilename + " saving...")
        np.savetxt(savefilename,
                   np.vstack((np.array(polyaIn) / 3, np.array(polyaIne) / 3,
                              inv * np.array(suspIn), inv * np.array(suspIne),
                              np.array(polyaOut) / 3, np.array(polyaOute) / 3,
                              outv * np.array(suspOut), outv * np.array(suspOute))), delimiter=",")
        return
        # return np.vstack((np.array(polyaIn) / 3, np.array(polyaIne) / 3,
        #                   inv * np.array(suspIn), inv * np.array(suspIne),
        #                   np.array(polyaOut) / 3, np.array(polyaOute) / 3,
        #                   outv * np.array(suspOut), outv * np.array(suspOute)))
    resvAll = []
    ressAll = []
    resvIn = []
    ressIn = []
    suspAll = []
    suspeAll = []
    suspIn = []
    suspeIn = []
    _, nums, sumnum, edgenum = getPrAndNumberOfSites(4 * nt + 2)
    for i in range(0, 11):
        arr = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
            diskName,
            getHeadByNt(nt),
            nt,
            folderNames[0],
            getHeadByNt(nt),
            getFileName(measureVal),
            nt,
            i))
        arr = arr * signs[0]
        if len(folderNames) > 1:
            for j in range(1, len(folderNames)):
                arrj = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
                    diskName,
                    getHeadByNt(nt),
                    nt,
                    folderNames[j],
                    getHeadByNt(nt),
                    getFileName(measureVal),
                    nt,
                    i))
                arrj = arrj * signs[j]
                arr = np.vstack((arr, arrj))
        rearAll = np.dot(np.real(arr), nums) / sumnum
        v, s = JacknifeMean(rearAll)
        _, _, t = AutoCorrelationSingleVariable(rearAll)
        resvAll.append(v)
        ressAll.append(s * np.sqrt(2 * t))
        if MeasureValueType.Chiral == measureVal:
            v, s = JacknifeCumulant(rearAll)
            suspAll.append(v)
            suspeAll.append(s * np.sqrt(2 * t))
        rearIn = np.dot(np.real(arr)[:,:edgenum], nums[:edgenum]) / np.sum(nums[:edgenum])
        v, s = JacknifeMean(rearIn)
        _, _, t = AutoCorrelationSingleVariable(rearIn)
        resvIn.append(v)
        ressIn.append(s * np.sqrt(2 * t))
        if MeasureValueType.Chiral == measureVal:
            v, s = JacknifeCumulant(rearIn)
            suspIn.append(v)
            suspeIn.append(s * np.sqrt(2 * t))
    if measureVal is MeasureValueType.Chiral:
        outv = getOutVolume(nt) * nt * nt
        inv = getInVolume(nt) * nt * nt
        factor = getKappa(nt) * 2
        savefilename = savefolder + getSaveName(measureVal) + "_nt{}.csv".format(nt)
        print(savefilename + " saving...")
        np.savetxt(savefilename, np.vstack((
            factor * np.array(resvIn), factor * np.array(ressIn),
            factor * factor * inv * np.array(suspIn), factor * factor * inv * np.array(suspeIn),
            factor * np.array(resvAll), factor * np.array(ressAll),
            factor * factor * outv * np.array(suspAll), factor * factor * outv * np.array(suspeAll))), delimiter=",")
        return
        # return np.vstack((
        #     factor * np.array(resvIn), factor * np.array(ressIn),
        #     factor * factor * inv * np.array(suspIn), factor * factor * inv * np.array(suspeIn),
        #     factor * np.array(resvAll), factor * np.array(ressAll),
        #     factor * factor * outv * np.array(suspAll), factor * factor * outv * np.array(suspeAll)))
    savefilename = savefolder + getSaveName(measureVal) + "_nt{}.csv".format(nt)
    print(savefilename + " saving...")
    if measureVal in [MeasureValueType.JFL, MeasureValueType.JFPot, MeasureValueType.JFS]:
        np.savetxt(savefilename,
                   np.vstack((np.array(resvIn) / 2, np.array(ressIn) / 2, np.array(resvAll) / 2, np.array(ressAll) / 2)),
                   delimiter=",")
        return
    np.savetxt(savefilename, np.vstack((np.array(resvIn), np.array(ressIn), np.array(resvAll), np.array(ressAll))), delimiter=",")
    # return np.vstack((np.array(resvIn), np.array(ressIn), np.array(resvAll), np.array(ressAll)))


def getRDist(diskName, folderNames, measureVal: MeasureValueType, nt: int, omegalst: list, savefolder: str):
    """
    Polyakov use real
    2 kappa is multiplied for chiral
    fermion momentum, Nf is divided
    """
    polyakovUseAbsOrReal = False
    resvall = []
    ressall = []
    pr, nums, _, _ = getPrAndNumberOfSites(4 * nt + 2)
    for omega in omegalst:
        if measureVal is not MeasureValueType.Polya and measureVal is not MeasureValueType.Chiral:
            if 0 == omega:
                continue
        resv = []
        ress = []
        arr = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
            diskName,
            getHeadByNt(nt),
            nt,
            folderNames[0],
            getHeadByNt(nt),
            getFileName(measureVal),
            nt,
            omega))
        if len(folderNames) > 1:
            for j in range(1, len(folderNames)):
                arrj = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
                    diskName,
                    getHeadByNt(nt),
                    nt,
                    folderNames[j],
                    getHeadByNt(nt),
                    getFileName(measureVal),
                    nt,
                    omega))
                arr = np.vstack((arr, arrj))
        if measureVal is MeasureValueType.Polya:
            rearAll = np.abs(arr) if polyakovUseAbsOrReal else np.real(arr)
            rearAll = rearAll / 3
        elif measureVal is MeasureValueType.Chiral:
            rearAll = np.real(arr) * 2 * getKappa(nt)
        elif measureVal in [MeasureValueType.JFPot, MeasureValueType.JFL, MeasureValueType.JFS]:
            rearAll = np.real(arr) / 2
        else:
            rearAll = np.real(arr)
        for i in range(len(pr)):
            if measureVal in [MeasureValueType.JG, MeasureValueType.JFL, MeasureValueType.JGChen, MeasureValueType.JGPot, MeasureValueType.JFPot, MeasureValueType.JGSurf]:
                if 0 == i:
                    resv.append(0)
                    ress.append(0)
                    continue
            v, s = JacknifeMean(rearAll[:, i])
            _, _, t = AutoCorrelationSingleVariable(rearAll[:, i])
            resv.append(v)
            ress.append(s * np.sqrt(2 * t))
        resvall.append(resv)
        ressall.append(ress)
    savefilename = savefolder + getSaveName(measureVal) + "_rdist_nt{}.csv".format(nt)
    print(savefilename + " saving...")
    np.savetxt(savefilename, np.vstack((np.array(pr), np.array(resvall), np.array(ressall))), delimiter=",")

def getRDistOneFile(pr, arr):
    resv = []
    ress = []
    for i in range(len(pr)):
        v, s = JacknifeMean(arr[:, i])
        _, _, t = AutoCorrelationSingleVariable(arr[:, i])
        resv.append(v)
        ress.append(s * np.sqrt(2 * t))
    return resv, ress

def getVDist(diskName, folderNames, measureVal: MeasureValueType, nt: int, edge: int):
    """
    For Polyakov, real is used
    2kappa is multiplied
    """
    polyakovUseAbsOrReal = False
    resx = []
    resv = []
    ress = []
    ressv = []
    resss = []
    pr, nums, _, edge = getPrAndNumberOfSitesTighter(4 * nt + 2, edge)
    for omega in range(0, 11):
        arr = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
            diskName,
            getHeadByNt(nt),
            nt,
            folderNames[0],
            getHeadByNt(nt),
            getFileName(measureVal),
            nt,
            omega))
        if len(folderNames) > 1:
            for j in range(1, len(folderNames)):
                arrj = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
                    diskName,
                    getHeadByNt(nt),
                    nt,
                    folderNames[j],
                    getHeadByNt(nt),
                    getFileName(measureVal),
                    nt,
                    omega))
                arr = np.vstack((arr, arrj))
        if measureVal is MeasureValueType.Polya:
            rearAll = np.abs(arr) if polyakovUseAbsOrReal else np.real(arr)
        else:
            rearAll = np.real(arr)
        for i in range(edge):
            if measureVal is MeasureValueType.Polya:
                v, s = JacknifeMean(rearAll[:, i])
                _, _, t = AutoCorrelationSingleVariable(rearAll[:, i])
                resx.append(omega * pr[i] * getOmegaSep(nt))
                resv.append(v / 3)
                ress.append(s * np.sqrt(2 * t) / 3)
                v, s = JacknifeCumulant(rearAll[:, i])
                # nums[i] is num[i] in x-y plane
                # (nums[i] * nt * 2) is spatial volume
                factor = nums[i] * nt * 2
                ressv.append(factor * v / 9)
                resss.append(factor * s * np.sqrt(2 * t) / 9)
            elif measureVal is MeasureValueType.Chiral:
                factorchiral = 2 * getKappa(nt)
                v, s = JacknifeMean(rearAll[:, i])
                _, _, t = AutoCorrelationSingleVariable(rearAll[:, i])
                resx.append(omega * pr[i] * getOmegaSep(nt))
                resv.append(factorchiral * v)
                ress.append(factorchiral * s * np.sqrt(2 * t))
                v, s = JacknifeCumulant(rearAll[:, i])
                # nums[i] is num[i] in x-y plane
                # (nums[i] * nt * 2) is spatial volume
                # (nums[i] * nt * 2) * nt is volume
                # (nums[i] * nt * 2) * nt * nt * nt is volume / T^2
                factor = nums[i] * nt * 2 * nt * nt * nt
                ressv.append(factorchiral * factorchiral * factor * v)
                resss.append(factorchiral * factorchiral * factor * s * np.sqrt(2 * t))
            else:
                print("should only support Polyakov and Chiral")
    return np.array(resx), np.array(resv), np.array(ress), np.array(ressv), np.array(resss)

def getRhoXi(diskName, folderNames, measureVal: MeasureValueType, nt: int):
    """
    xi is considered with r=0 to edge
    rho is considered with r=1 to edge
    omega start from 1
    For fermion angular momenta, Nf is divided.
    """
    resv = []
    ress = []
    pr, nums, _, edge = getPrAndNumberOfSites(4 * nt + 2)
    volume = np.sum(nums[1:edge])
    pr = np.array(pr)[1:edge]
    prSq = pr * pr
    for omega in range(1, 11):
        arr = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
            diskName,
            getHeadByNt(nt),
            nt,
            folderNames[0],
            getHeadByNt(nt),
            getFileName(measureVal),
            nt,
            omega))
        if len(folderNames) > 1:
            for j in range(1, len(folderNames)):
                arrj = LoadMathematicaCSV("{}{}Nt{}\\{}\\NC{}__{}_Nt{}_O{}.csv".format(
                    diskName,
                    getHeadByNt(nt),
                    nt,
                    folderNames[j],
                    getHeadByNt(nt),
                    getFileName(measureVal),
                    nt,
                    omega))
                arr = np.vstack((arr, arrj))
        rearAll = np.real(arr)
        if measureVal is MeasureValueType.JFS or measureVal is MeasureValueType.JFL:
            rearAll = rearAll / 2 # remove Nf
        if measureVal is MeasureValueType.JGS or measureVal is MeasureValueType.JFS:
            rearAll = rearAll[:, 0:edge]
            rearAll = np.dot(rearAll, nums[0:edge]) / (omega * getOmegaSep(nt)) / volume
        else:
            rearAll = rearAll[:, 1:edge]
            rearAll = rearAll / prSq
            rearAll = np.dot(rearAll, nums[1:edge]) / (omega * getOmegaSep(nt)) / volume
        v, s = JacknifeMean(rearAll)
        _, _, t = AutoCorrelationSingleVariable(rearAll)
        resv.append(v)
        ress.append(s * np.sqrt(2 * t))
    return np.array(resv), np.array(ress)
