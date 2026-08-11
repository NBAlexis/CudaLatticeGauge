import numpy as np

from JacknifePrograms import LoadMathematicaCSV
from JointCorrelatedFit import JointCorrelatedFitAndDraw

data1 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUp\\TUp6__565_bosonvalueR_ZSlice.csv"))
data2 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TDown\\TDown6__565_bosonvalueR_ZSlice.csv"))
data3 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUD\\TUD6__565_bosonvalueR_ZSlice.csv"))
data4 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUDA\\TUDA6__565_bosonvalueR_ZSlice.csv"))
data5 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUp\\TUp6__615_bosonvalueR_ZSlice.csv"))
data6 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TDown\\TDown6__615_bosonvalueR_ZSlice.csv"))
data7 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUD\\TUD6__615_bosonvalueR_ZSlice.csv"))
data8 = np.real(LoadMathematicaCSV("I:\\BetaGradient\\TemperatureDist\\TUDA\\TUDA6__615_bosonvalueR_ZSlice.csv"))

xdata1 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata2 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata3 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata4 = np.array([n + 1 for n in range(7)] + [n + 9 for n in range(15)])
xdata5 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata6 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata7 = np.array([n + 1 for n in range(11)] + [n + 13 for n in range(11)])
xdata8 = np.array([n + 1 for n in range(7)] + [n + 9 for n in range(15)])

ydata1 = data1[:, xdata1]
ydata2 = data2[:, xdata2]
ydata3 = data3[:, xdata3]
ydata4 = data4[:, xdata4]
ydata5 = data5[:, xdata5]
ydata6 = data6[:, xdata6]
ydata7 = data7[:, xdata7]
ydata8 = data8[:, xdata8]

initialguess = [1.0]
bounds = [[0.5, 1.5]]
# ['L-BFGS-B', 'TNC', 'SLSQP']
methods = ""

def func1(param, x):
    return 1.14538181 + (2.29076362 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 12)))

def func2(param, x):
    return 1.14538181 + (0.5726909 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 12)))

def func3(param, x):
    return 1.14538181 + (2.0 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 12))) + (0.5 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x))) + (0.5 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 24)))

def func4(param, x):
    return 1.14538181 + (2.0 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 8))) + (0.5 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x))) + (0.5 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 24)))

def func5(param, x):
    return 1.14378345 + (2.2875669 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 12)))

def func6(param, x):
    return 1.14378345 + (0.57189173 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 12)))

def func7(param, x):
    return 1.14378345 + (2.0 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 12))) + (0.5 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x))) + (0.5 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 24)))

def func8(param, x):
    return 1.14378345 + (2.0 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 8))) + (0.5 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x))) + (0.5 - 1.14378345) * (1-np.tanh(param[0] * np.abs(x - 24)))

JointCorrelatedFitAndDraw([xdata1, xdata2, xdata3, xdata4, xdata5, xdata6, xdata7, xdata8],
                          [ydata1, ydata2, ydata3, ydata4, ydata5, ydata6, ydata7, ydata8],
                           [func1, func2, func3, func4, func5, func6, func7, func8], initialguess, bounds, methods)

JointCorrelatedFitAndDraw([xdata1, xdata2, xdata3, xdata4],
                          [ydata1, ydata2, ydata3, ydata4],
                           [func1, func2, func3, func4], initialguess, bounds, methods)

JointCorrelatedFitAndDraw([xdata5, xdata6, xdata7, xdata8],
                          [ydata5, ydata6, ydata7, ydata8],
                           [func5, func6, func7, func8], initialguess, bounds, methods)


"""

============================================================
JOINT CORRELATED FIT RESULTS
============================================================
Success: True
Message: Optimization terminated successfully.

Fit Quality Statistics:
  Total χ²: 342.7547
  Degrees of freedom: 176
  χ²/dof: 1.9475
  Data points: 176
  Parameters: 0
  Quality: ACCEPTABLE

χ² per dataset:
  Dataset 0: χ² = 98.9543, dof = 22, χ²/dof = 4.4979
  Dataset 1: χ² = 11.1570, dof = 22, χ²/dof = 0.5071
  Dataset 2: χ² = 60.1684, dof = 22, χ²/dof = 2.7349
  Dataset 3: χ² = 60.3586, dof = 22, χ²/dof = 2.7436
  Dataset 4: χ² = 29.1247, dof = 22, χ²/dof = 1.3238
  Dataset 5: χ² = 14.5125, dof = 22, χ²/dof = 0.6597
  Dataset 6: χ² = 34.3883, dof = 22, χ²/dof = 1.5631
  Dataset 7: χ² = 34.0910, dof = 22, χ²/dof = 1.5496

Fitted parameters with errors:
----------------------------------------------------------------------
Param_0    =   1.092090 ± 0.004764

============================================================
JOINT CORRELATED FIT RESULTS
============================================================
Success: True
Message: Optimization terminated successfully.

Fit Quality Statistics:
  Total χ²: 197.8295
  Degrees of freedom: 88
  χ²/dof: 2.2481
  Data points: 88
  Parameters: 0
  Quality: POOR

χ² per dataset:
  Dataset 0: χ² = 64.0180, dof = 22, χ²/dof = 2.9099
  Dataset 1: χ² = 21.3369, dof = 22, χ²/dof = 0.9699
  Dataset 2: χ² = 56.4279, dof = 22, χ²/dof = 2.5649
  Dataset 3: χ² = 56.0467, dof = 22, χ²/dof = 2.5476

Fitted parameters with errors:
----------------------------------------------------------------------
Param_0    =   1.053809 ± 0.006479

============================================================
JOINT CORRELATED FIT RESULTS
============================================================
Success: True
Message: Optimization terminated successfully.

Fit Quality Statistics:
  Total χ²: 81.8723
  Degrees of freedom: 88
  χ²/dof: 0.9304
  Data points: 88
  Parameters: 0
  Quality: GOOD

χ² per dataset:
  Dataset 0: χ² = 25.6184, dof = 22, χ²/dof = 1.1645
  Dataset 1: χ² = 8.0544, dof = 22, χ²/dof = 0.3661
  Dataset 2: χ² = 24.2877, dof = 22, χ²/dof = 1.1040
  Dataset 3: χ² = 23.9119, dof = 22, χ²/dof = 1.0869

Fitted parameters with errors:
----------------------------------------------------------------------
Param_0    =   1.129492 ± 0.007008

"""