import numpy as np

from AutoCorrelation import AutoCorrelation, AutoCorrelationSingleVariable

lst1 = 1.0 + np.random.uniform(-0.3, 0.3, 10000)
lst2 = 2.0 + np.random.uniform(-0.1, 0.1, 10000)
datatest = np.transpose(np.array([lst1, lst2]))
print(np.shape(datatest))

def functest(lst):
    return lst[0] + lst[1] * lst[1]


print(AutoCorrelation(datatest, functest, debugInfo=True))

lst3 = 3.0 + np.random.uniform(-0.3, 0.3, 10000)
print(AutoCorrelationSingleVariable(lst3, debugInfo=True))