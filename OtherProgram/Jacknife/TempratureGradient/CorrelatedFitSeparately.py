import numpy as np

from CorrelatedFit import CorrelatedFitAndDraw
from JacknifePrograms import LoadMathematicaCSV

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
methods = "SLSQP"

"""
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

"""

def func1(param, x):
    return 1.14538181 + (2.29076362 - 1.14538181) * (1-np.tanh(param[0] * np.abs(x - 12)))

CorrelatedFitAndDraw(xdata1, ydata1, func1, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata2, ydata2, func2, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata3, ydata3, func3, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata4, ydata4, func4, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata5, ydata5, func5, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata6, ydata6, func6, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata7, ydata7, func7, initialguess, bounds, methods)
CorrelatedFitAndDraw(xdata8, ydata8, func8, initialguess, bounds, methods)

"""
=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 1.8870
Fitted parameters: [0.99071446], Parameter errors (Bootstrap): [0.00080502]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 0.4430
Fitted parameters: [1.11814257], Parameter errors (Bootstrap): [0.00068319]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 2.6375
Fitted parameters: [1.06573408], Parameter errors (Bootstrap): [0.00109579]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 2.6254
Fitted parameters: [1.06485128], Parameter errors (Bootstrap): [0.00082171]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 1.1955
Fitted parameters: [1.11943463], Parameter errors (Bootstrap): [0.00070154]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 0.3610
Fitted parameters: [1.14329107], Parameter errors (Bootstrap): [0.00060284]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 1.1561
Fitted parameters: [1.13068966], Parameter errors (Bootstrap): [0.00071678]
fit using : SLSQP

=== Fit Results ===
Success: True
Message: Optimization terminated successfully
Chi2/dof: 1.1382
Fitted parameters: [1.13069277], Parameter errors (Bootstrap): [0.00066462]

"""