import numpy as np

"""
Implementation of 10.1016/s0010-4655(03)00467-3 with R=1
"""


def correlationAB(data, mean, a, b, t):
    """
    Equation (31)
    """
    l = len(data) - t
    data1 = data[0:l, a]
    data2 = data[t:(l + t), b]
    return np.dot(data1 - mean[a], data2 - mean[b]) / l


def fa(data, mean, a, func):
    """
    Equation (39)
    ha is defined in Equation (38)
    """
    ha = np.sqrt(correlationAB(data, mean, a, a, 0) / len(data))
    left = mean.copy()
    right = mean.copy()
    left[a] = left[a] + ha
    right[a] = right[a] - ha
    return (func(left) - func(right)) / (2 * ha)


def combineVariables(fa, data):
    """
    Equation (37)
    """
    return np.dot(data, fa)


def cw(cf, W, n):
    """
    Equation (35) and (49)
    """
    return (cf[0] if 0 == W else (cf[0] + 2 * np.sum(cf[1:W]))) * (1 + (2 * W + 1) / n)


def tau(W, cw, s):
    """
    Equation (51)

    where in (51), tau_int_F is in (41): tau_int_F = cF(W) / 2vF
    where vF is in (34):  vF = GammaF(0) (here we use the enlarged cF(0))
    so, (2tau_int_F+1)/(2tau_int_F-1) = (2cF(W)+2vF)/(2cF(W)-2vF) = (cF(W)+cF(0))/(cF(W)-vF(0))
    """
    cwsum = cw[W] + cw[0]
    cwdif = max(1.0e-12, cw[W] - cw[0])
    return s / np.log(cwsum / cwdif)


def gw(W, cw, s, n):
    """
    Equation (52)
    """
    tauW = tau(W, cw, s)
    return np.exp(-W / tauW) - tauW / np.sqrt(W * n)


def AutoCorrelation(data, func, s=1.5, debugInfo=False):
    dataCount = len(data)
    mean = np.mean(data, axis=0)
    combdata = combineVariables(np.array([fa(data, mean, i, func) for i in range(0, len(mean))]), data)
    combmean = np.array([np.mean(combdata)])
    combdata = np.transpose(np.array([combdata]))
    """
    starting from t=0
    """
    t = 0
    cflst = np.array([correlationAB(combdata, combmean, 0, 0, t)])
    cwlst = np.array([cw(cflst, t, dataCount)])
    lastGW = 1
    while lastGW > 0:
        t = t + 1
        cflst = np.append(cflst, correlationAB(combdata, combmean, 0, 0, t))
        cwlst = np.append(cwlst, cw(cflst, t, dataCount))
        lastGW = gw(t, cwlst, s, dataCount)
        if debugInfo:
            print("W: %d, G: %f" % (t, lastGW))
    finalCW = cwlst[t]
    """
    error is in Equation (44),
    tau independent is in Equation (41)
    """
    return func(mean), np.sqrt(finalCW / dataCount), max(0.5, finalCW / (2 * cwlst[0]))


"""
This is for the most cases where data is single variable
"""


def correlationABSingleVariable(data, mean, t):
    l = len(data) - t
    data1 = data[0:l]
    data2 = data[t:(l + t)]
    return np.dot(data1 - mean, data2 - mean) / l


def AutoCorrelationSingleVariable(data, s=1.5, debugInfo=False):
    dataCount = len(data)
    mean = np.mean(data)
    t = 0
    cflst = np.array([correlationABSingleVariable(data, mean, t)])
    cwlst = np.array([cw(cflst, t, dataCount)])
    lastGW = 1
    while lastGW > 0:
        t = t + 1
        cflst = np.append(cflst, correlationABSingleVariable(data, mean, t))
        cwlst = np.append(cwlst, cw(cflst, t, dataCount))
        lastGW = gw(t, cwlst, s, dataCount)
        if debugInfo:
            print("W: %d, G: %f" % (t, lastGW))
    finalCW = cwlst[t]
    return mean, np.sqrt(finalCW / dataCount), max(0.5, finalCW / (2 * cwlst[0]))
