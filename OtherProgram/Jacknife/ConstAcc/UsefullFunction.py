from enum import Enum

import numpy as np

from AutoCorrelation import AutoCorrelationSingleVariable
from JacknifePrograms import LoadMathematicaCSV, JacknifeMean, JacknifeCumulant, PrintAsMatlabArray


class MeasureValueType(Enum):
    POLYA = 0
    CHIRAL = 1
    G3 = 2
    G4 = 3

def ReCalculateCondensate(folderName, saveFolder, mid, measureType, startIdx):
    saveName = ""
    if MeasureValueType.CHIRAL == measureType:
        saveName = "chiral"
    elif MeasureValueType.G3 == measureType:
        saveName = "g3"
    elif MeasureValueType.G4 == measureType:
        saveName = "g4"
    maxz = 12
    lst = ["m40", "m35", "m30", "m25", "m20", "m15", "m10", "m05", "00", "05", "10", "15", "20", "25", "30", "35", "40"]
    glst = [-0.04 + 0.005 * i for i in range(17)]
    if mid:
        saveName = "mid" + saveName
        maxz = 13
        lst = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
        glst = [0 + 0.01 * i for i in range(11)]
    lst3 = ["56", "565", "57", "575", "58", "585", "59", "595", "60", "605"]
    allchiral = []
    allchirale = []
    for headname in lst3:
        print("dealing with " + headname)
        chiral = []
        chirale = []
        for j in range(len(lst)):
            filename = folderName.format(headname, headname, lst[j])
            arr = LoadMathematicaCSV(filename)[startIdx:]
            if MeasureValueType.CHIRAL == measureType:
                arr = -np.real(arr)
            elif MeasureValueType.G3 == measureType:
                arr = -np.imag(arr)
            elif MeasureValueType.G4 == measureType:
                arr = np.imag(arr)
            g = glst[j]
            denorm = np.array([1 + g * z for z in range(0, maxz)])
            if mid:
                denorm = np.array([1 + g * (z - 6) for z in range(0, maxz)])
                arr = arr * np.array([1 + g * (z - 6) for z in range(0, maxz)])
            else:
                arr = arr * np.array([1 + g * z for z in range(0, maxz)])
            arr = np.sum(arr, axis=1) / np.sum(denorm)
            v, s = JacknifeMean(arr)
            _, _, t = AutoCorrelationSingleVariable(arr)
            chiral.append(v)
            chirale.append(s * np.sqrt(2 * t))
        allchiral.append(chiral)
        allchirale.append(chirale)
    allchiral = np.array(allchiral)
    allchirale = np.array(allchirale)
    np.savetxt(saveFolder + saveName + ".csv", np.vstack((allchiral, allchirale)), delimiter=',')
    print(saveFolder + saveName + ".csv" + " saved successfully!")

def Export(folderName, saveFolder, mid, measureType, startidx):
    if measureType != MeasureValueType.POLYA:
        ReCalculateCondensate(folderName, saveFolder, mid, measureType, startidx)
        return
    headnames = ["56", "565", "57", "575", "58", "585", "59", "595", "60", "605"]
    lst = ["m40", "m35", "m30", "m25", "m20", "m15", "m10", "m05", "00", "05", "10", "15", "20", "25", "30", "35", "40"]
    saveName = "polya"
    if mid:
        saveName = "mid" + saveName
        lst = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
    alltosave = None
    for headname in headnames:
        polyalst = []
        polyalste = []
        susplst = []
        susplste = []
        times = []
        for j in range(len(lst)):
            filename = folderName.format(headname, headname, lst[j])
            arr = LoadMathematicaCSV(filename)[startidx:]
            arr = np.abs(arr)
            v, s = JacknifeMean(arr, headname + lst[j])
            polyalst.append(v)
            polyalste.append(s)
            v, s = JacknifeCumulant(arr)
            susplst.append(v)
            susplste.append(s)
            _, _, t = AutoCorrelationSingleVariable(arr)
            times.append(t)
        polyalst = np.array(polyalst)
        polyalste = np.array(polyalste)
        susplst = np.array(susplst)
        susplste = np.array(susplste)
        times = np.array(times)
        polyalste = polyalste * np.sqrt(2 * times)
        if alltosave is None:
            alltosave = polyalst
        else:
            alltosave = np.vstack((alltosave, polyalst))
        alltosave = np.vstack((alltosave, polyalste))
        alltosave = np.vstack((alltosave, susplst))
        susplste = susplste * np.sqrt(2 * times)
        alltosave = np.vstack((alltosave, susplste))
    np.savetxt(saveFolder + saveName + ".csv", alltosave, delimiter=',')
    print(saveFolder + saveName + ".csv" + " saved successfully!")

def ExportZSlice(folderName, saveFolder, mid, measureType, startidx):
    maxZ = 12
    lst = ["m40", "m35", "m30", "m25", "m20", "m15", "m10", "m05", "00", "05", "10", "15", "20", "25", "30", "35", "40"]
    lst3 = ["56", "565", "57", "575", "58", "585", "59", "595", "60", "605"]
    saveName = "polya"
    if MeasureValueType.CHIRAL == measureType:
        saveName = "chiral"
    elif MeasureValueType.G3 == measureType:
        saveName = "g3"
    elif MeasureValueType.G4 == measureType:
        saveName = "g4"
    if mid:
        saveName = "mid" + saveName
        maxZ = 13
        lst = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
    for headname in lst3:
        ptest = []
        petest = []
        susp = []
        suspe = []
        tset = []
        print("dealing with " + headname)
        for j in range(len(lst)):
            filename = folderName.format(headname, headname, lst[j])
            arr = LoadMathematicaCSV(filename)[startidx:]
            if MeasureValueType.POLYA == measureType:
                arr = np.abs(arr)
            elif MeasureValueType.CHIRAL == measureType:
                """
                arr was -Tr[D]
                we need Tr[D] which is -real(arr)
                """
                arr = -np.real(arr)
            elif MeasureValueType.G3 == measureType:
                """
                arr was -Tr[D]
                we need i * (-Tr[D]) which is -i Tr[D] which is -imag(arr)
                """
                arr = -np.imag(arr)
            elif MeasureValueType.G4 == measureType:
                """
                arr was -Tr[D]
                we need Im[-Tr[D]] which is i * Tr[D] which is imag(arr)
                """
                arr = np.imag(arr)
            polyalst = []
            polyalste = []
            susplst = []
            susplste = []
            times = []
            for z in range(0, maxZ):
                thisZSlice = arr[:, z]
                thisZSlice = thisZSlice.flatten()
                v, s = JacknifeMean(thisZSlice, headname + lst[j] if z == 0 else "")
                polyalst.append(v)
                polyalste.append(s)
                if MeasureValueType.POLYA == measureType:
                    v, s = JacknifeCumulant(thisZSlice)
                    susplst.append(v)
                    susplste.append(s)
                v, s, t = AutoCorrelationSingleVariable(thisZSlice)
                times.append(t)
            ptest.append(polyalst)
            petest.append(polyalste)
            if MeasureValueType.POLYA == measureType:
                susp.append(susplst)
                suspe.append(susplste)
            tset.append(times)
        ptest = np.array(ptest)
        petest = np.array(petest)
        if MeasureValueType.POLYA == measureType:
            susp = np.array(susp)
            suspe = np.array(suspe)
        tset = np.array(tset)
        petest = petest * np.sqrt(2 * tset)
        if MeasureValueType.POLYA == measureType:
            suspe = suspe * np.sqrt(2 * tset)
        np.savetxt(saveFolder + saveName + headname + ".csv", ptest, delimiter=',')
        np.savetxt(saveFolder + saveName + "e" + headname + ".csv", petest, delimiter=',')
        if MeasureValueType.POLYA == measureType:
            np.savetxt(saveFolder + saveName + "s" + headname + ".csv", susp, delimiter=',')
            np.savetxt(saveFolder + saveName + "se" + headname + ".csv", suspe, delimiter=',')
