import lmfit
from lmfit import Model, Parameter, report_fit
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, fsolve

from MesonStructures import all_signs
from UsefulFunctions import jackknife_statistic

type = 16
all_signs_lst = all_signs()
all_type = len(all_signs_lst[type])

correlation_func = np.load(f"data/correlationpp_{type}_0.npy")
for i in range(1, all_type):
    correlation_func = correlation_func + np.load(f"data/correlationpp_{type}_{i}.npy")
# print(correlation_func.tolist())
print(np.shape(correlation_func))


def solve_params(x_data, y_data, L, initial_guess=None):
    x1, x2, x3, x4 = map(float, x_data)
    y1, y2, y3, y4 = map(float, y_data)

    if initial_guess is None:
        y_avg = np.mean(y_data)
        y_range = np.max(y_data) - np.min(y_data)
        A1_guess = y_avg / 2
        A2_guess = y_range / 4
        m1_guess = 0.1
        m2_guess = 0.1
        initial_guess = [A1_guess, A2_guess, m1_guess, m2_guess]

    def equations(params):
        A1, A2, m1, m2 = params

        eq1 = A1 * (np.exp(-m1 * x1) + np.exp(-m1 * (L - x1))) + A2 * np.cos(np.pi * x1) * (
                    np.exp(-m2 * x1) + np.exp(-m2 * (L - x1))) - y1
        eq2 = A1 * (np.exp(-m1 * x2) + np.exp(-m1 * (L - x2))) + A2 * np.cos(np.pi * x2) * (
                    np.exp(-m2 * x2) + np.exp(-m2 * (L - x2))) - y2
        eq3 = A1 * (np.exp(-m1 * x3) + np.exp(-m1 * (L - x3))) + A2 * np.cos(np.pi * x3) * (
                    np.exp(-m2 * x3) + np.exp(-m2 * (L - x3))) - y3
        eq4 = A1 * (np.exp(-m1 * x4) + np.exp(-m1 * (L - x4))) + A2 * np.cos(np.pi * x4) * (
                    np.exp(-m2 * x4) + np.exp(-m2 * (L - x4))) - y4
        return [eq1, eq2, eq3, eq4]
    solution, info, ier, msg = fsolve(equations, initial_guess, full_output=True)
    if ier != 1:
        print(f"Warning, solver did not converge：{msg}")
    return tuple(solution)

real_res = jackknife_statistic(np.real(correlation_func))

# plt.plot(range(48), real_res['std_error'])
# plt.plot(range(48), np.abs(np.mean(np.imag(correlation_func), axis=0)))
# plt.show()

x_data = np.array([t for t in range(48)])
y_data = real_res['jackknife_estimate']

x_data2 = []
y_data2 = []
for t2 in range(1, 24):
    x_data2.append(t2*1.0)
    y_data2.append((y_data[t2] + y_data[48 - t2])/2)

"""
for tstart in range(1, 12):
    x_data_to_fit = x_data2[tstart:(tstart + 4)]
    y_data_to_fit = y_data2[tstart:(tstart + 4)]
    print(solve_params(x_data_to_fit, y_data_to_fit, 48))
"""

def my_func(t, Ap, mp, Am, mm):
    return Ap * (np.exp(-mp * t) + np.exp(-mp * (48 - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (48 - t)))

def my_func2(t, Ap, mp, Am, mm):
    return Ap * (np.exp(-mp * t) + np.exp(-mp * (48 - t))) + np.cos(np.pi * t) * Am * (np.exp(-mm * t) + np.exp(-mm * (48 - t)))

model = Model(my_func)
params = lmfit.Parameters()
params.add('Ap', value=np.mean(y_data2), vary=True, min=-np.inf, max=np.inf)
params.add('mp', value=0.1, vary=True, min=0.00001, max=10.0)
params.add('Am', value=np.std(y_data2), vary=True, min=-np.inf, max=np.inf)
params.add('mm', value=0.2, vary=True, min=0.00001, max=10.0)

tsep = 8
for tstart in range(8):
    result = model.fit(y_data2[tstart:(tstart+tsep)], params, t=x_data2[tstart:(tstart+tsep)])

    print(result.fit_report())
    params_fitted = [result.params['Ap'].value,
                     result.params['mp'].value,
                     result.params['Am'].value,
                     result.params['mm'].value]


"""
bounds = ([-np.inf, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf])
params = curve_fit(my_func, x_data2[fit_start_t:fit_end_t], y_data2[fit_start_t:fit_end_t], maxfev=50000, bounds=bounds)
print("res:", params[0])
"""

# """
xx_data = []
yy_data = []
for t in range(48):
    xx_data.append(t)
    yy_data.append(my_func(t, result.params['Ap'].value, result.params['mp'].value, result.params['Am'].value, result.params['mm'].value))

plt.plot(x_data[1:], y_data[1:])
plt.plot(x_data2, y_data2)
plt.plot(xx_data[1:], yy_data[1:], 'o')
plt.show()
# """