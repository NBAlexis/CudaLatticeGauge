from io import StringIO

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def LoadMathematicaCSV(fileName: str):
    with open(fileName) as f:
        contents = f.read()
        contents = contents.replace("*^", "e")
        contents = contents.replace(" I", "j")
        contents = contents.replace(" ", "")
        c = StringIO(contents)
        d = np.loadtxt(c, delimiter=",", dtype=str)
        try:
            d = d.astype(complex)
        except:
            for i in range(len(d)):
                try:
                    ds = d[i].astype(complex)
                except:
                    print(d[i])
        d = d.astype(complex)
        return d


def PrintProgressBar(title, d, n):
    prog = 100 * d / n
    i = int(prog)
    num = i // 2
    if 100 == prog:
        process = "\r%s: %d/%d [%3s%%]: |%-50s|\n" % (title, d, n, i, '|' * num)
    else:
        process = "\r%s: %d/%d [%3s%%]: |%-50s|" % (title, d, n, i, '|' * num)
    print(process, end='', flush=True)


def JacknifeMean(lst, progressBarTitle=""):
    size = len(lst)
    orignalv = np.mean(lst)
    argsv = 0
    avers = []
    for i in range(size):
        subseto = np.delete(lst, i)
        avers.append(np.mean(subseto))
        argsv = argsv + (avers[i] - orignalv) ** 2
        if 0 != len(progressBarTitle):
            PrintProgressBar(progressBarTitle, i, size)
    argsv = argsv * (size - 1) / size
    unbaisv = orignalv + (size - 1) * (orignalv - np.mean(np.array(avers)))
    if 0 != len(progressBarTitle):
        PrintProgressBar(progressBarTitle, size, size)
    return unbaisv, np.sqrt(argsv)


def JacknifeCumulant(lst, progressBarTitle=""):
    size = len(lst)
    orignalv = np.mean(lst * lst) - np.mean(lst) * np.mean(lst)
    argsv = 0
    avers = []
    for i in range(size):
        subseto = np.delete(lst, i)
        avers.append(np.mean(subseto * subseto) - np.mean(subseto) * np.mean(subseto))
        argsv = argsv + (avers[i] - orignalv) ** 2
        if 0 != len(progressBarTitle):
            PrintProgressBar(progressBarTitle, i, size)
    argsv = argsv * (size - 1) / size
    unbaisv = orignalv + (size - 1) * (orignalv - np.mean(np.array(avers)))
    if 0 != len(progressBarTitle):
        PrintProgressBar(progressBarTitle, size, size)
    return unbaisv, np.sqrt(argsv)


def expfit(x, cr, vr):
    return cr * np.exp(-vr * x)


def potentialfit(x, a, b, c):
    return a + b * x - c / x


def FitWilsonLoop(vr, pr, halfT, tStart, tEnd, maxR, showFig=False):
    pr = np.abs(pr)
    meanAbsR = np.abs(np.mean(vr, axis=0))
    meanWrt = [[meanAbsR[r * halfT + t] for t in range(halfT)] for r in range(len(pr))]
    crlst = []
    vrlst = []
    for r in range(len(pr)):
        xdata = [t + 1 for t in range(tStart, tEnd)]
        ydata = [meanWrt[r][t] for t in range(tStart, tEnd)]
        parameters, _ = curve_fit(expfit, xdata, ydata)
        crlst.append(parameters[0])
        vrlst.append(parameters[1])
    vrlst = np.array(vrlst)
    parameters2, _ = curve_fit(potentialfit, pr[0:maxR], vrlst[0:maxR])
    if showFig:
        fig = plt.figure()
        ax1 = fig.subplots()
        ax2 = ax1.twinx()
        rtoplot = np.array([0.05 + 0.1 * i for i in range(2 + int(10 * max(pr)))])
        ax1.plot(pr, vrlst, '+')
        ax1.plot(rtoplot, potentialfit(rtoplot, parameters2[0], parameters2[1], parameters2[2]))
        ax1.set_ylim((0, max(vrlst) + 0.1))
        ax2.plot(pr, crlst, 'o')
        ax1.set_xlabel('r[a]')
        ax1.set_ylabel('v(r)')
        ax2.set_ylabel('c(r)')
        plt.show()
    return parameters2[0], parameters2[1], parameters2[2]


def JacknifeWilsonLoop(vr, pr, halfT, tStart, tEnd, maxR, showProgress="", showFirtFig=False):
    aOrignal, bOrignal, cOrignal = FitWilsonLoop(vr, pr, halfT, tStart, tEnd, maxR, showFirtFig)
    orignalR0 = np.sqrt((1.65 - cOrignal) / bOrignal)
    orignalR1 = np.sqrt((1 - cOrignal) / bOrignal)
    orignalC0 = - aOrignal + cOrignal / (1.5 * orignalR0) - np.pi / (12 * 1.5 * orignalR0)
    argsr0 = 0
    argsr1 = 0
    argsc0 = 0
    meanr0 = 0
    meanr1 = 0
    meanc0 = 0
    size = len(vr)
    for i in range(size):
        vrreplace = np.delete(vr, i, 0)
        a, b, c = FitWilsonLoop(vrreplace, pr, halfT, tStart, tEnd, maxR)
        newR0 = np.sqrt((1.65 - c) / b)
        newR1 = np.sqrt((1 - c) / b)
        newC0 = - a + c / (1.5 * newR0) - np.pi / (12 * 1.5 * newR0)
        argsr0 = argsr0 + (newR0 - orignalR0) ** 2
        argsr1 = argsr1 + (newR1 - orignalR1) ** 2
        argsc0 = argsc0 + (newC0 - orignalC0) ** 2
        meanr0 = meanr0 + newR0
        meanr1 = meanr1 + newR1
        meanc0 = meanc0 + newC0
        if len(showProgress) > 0:
            PrintProgressBar(showProgress, i, size)
    argsr0 = np.sqrt((size - 1) * argsr0 / size)
    argsr1 = np.sqrt((size - 1) * argsr1 / size)
    argsc0 = np.sqrt((size - 1) * argsc0 / size)
    unbaisr0 = orignalR0 + (size - 1) * (orignalR0 - meanr0 / size)
    unbaisr1 = orignalR1 + (size - 1) * (orignalR1 - meanr1 / size)
    unbaisc0 = orignalC0 + (size - 1) * (orignalC0 - meanc0 / size)
    if len(showProgress) > 0:
        PrintProgressBar(showProgress, size, size)
    return unbaisr0, argsr0, unbaisr1, argsr1, unbaisc0, argsc0


def PrintAsMathematicaArray(arr, header="") -> str:
    ret = str(np.array(arr))
    ret = ret.replace("\n", "")
    ret = ret.replace("\r", "")
    ret = ret.replace(" ", ",")
    ret = ret.replace("j", "I")
    oldLen = len(ret)
    ret = ret.replace(",,", ",")
    while len(ret) != oldLen:
        oldLen = len(ret)
        ret = ret.replace(",,", ",")
    ret = ret.replace("[", "{")
    ret = ret.replace("]", "}")
    ret = ret.replace(",}", "}")
    ret = ret.replace("{,", "{")
    ret = ret.replace("e", "*^")
    ret = ret + ";"
    if 0 != len(header):
        ret = header + "=" + ret
    return ret


def PrintAsMatlabArray(arr, header="") -> str:
    ret = str(np.array(arr))
    ret = ret.replace("\n", "")
    ret = ret.replace("\r", "")
    ret = ret.replace(" ", ",")
    ret = ret.replace("j", "i")
    oldLen = len(ret)
    ret = ret.replace(",,", ",")
    while len(ret) != oldLen:
        oldLen = len(ret)
        ret = ret.replace(",,", ",")
        ret = ret.replace(" i", "i")
    ret = ret.replace(",]", "]")
    ret = ret.replace("[,", "[")
    ret = ret + ";"
    if 0 != len(header):
        ret = header + "=" + ret
    return ret
