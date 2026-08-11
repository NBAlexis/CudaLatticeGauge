import lmfit
from lmfit import Model, Parameter, report_fit
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, fsolve

from CorrelatedFit import correlated_fit_with_flexible_bounds
from MesonStructures import all_signs

nt = 48
t_stride = 8
t_end = 0
type = 0
all_signs_lst = all_signs()
all_type = len(all_signs_lst[type])

def my_func(t, args):
    Ap = args[0]
    mp = args[1]
    Am = args[2]
    mm = args[3]
    return Ap * (np.exp(-mp * t) + np.exp(-mp * (48 - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (48 - t)))

correlation_func = np.load(f"data/correlationpp_{type}_0.npy")
for i in range(1, all_type):
    correlation_func = correlation_func + np.load(f"data/correlationpp_{type}_{i}.npy")

correlation_func = np.real(correlation_func)
# pick a separation of t
final_res1 = []
final_res1e = []
final_res2 = []
final_res2e = []
check_step = (nt // 2) - t_stride - t_end
for t_start in range(1, check_step):
    t_value = np.array([k for k in range(nt)])
    pickidx = [k + t_start + 1 for k in range(t_stride)] + [nt - t_stride - t_start + k for k in range(t_stride)]
    x_pick = t_value[pickidx]
    y_pick = correlation_func[:, pickidx]
    initial_guess = np.array([np.mean(correlation_func), 0.1, np.std(correlation_func), 0.1])
    upper_bounds = [np.inf, 2.0, np.inf, 2.0]
    lower_bounds = [-np.inf, 0.1, -np.inf, 0.1]

    fitres = correlated_fit_with_flexible_bounds(
        y_pick,
        x_pick,
        my_func,
        initial_guess,
        (lower_bounds, upper_bounds))
    print(fitres['params'])
    print(fitres['params_err'])
    print(fitres['chi2'])
    print(fitres['chi2_dof'])
    final_res1.append(fitres['params'][1])
    final_res1e.append(fitres['params_err'][1])
    final_res2.append(fitres['params'][3])
    final_res2e.append(fitres['params_err'][3])
    # print(y_pick)

print(final_res1)
print(final_res1e)
print(final_res2)
print(final_res2e)

plt.errorbar([k + 1 for k in range(check_step - 1)], final_res2, final_res2e,fmt=':o')
plt.xlabel("$t_s$")
plt.ylabel("$am_{\pi}$")
plt.show()

plt.errorbar([k + 1 for k in range(check_step - 1)], final_res1, final_res1e,fmt=':o')
plt.xlabel("$t_s$")
plt.ylabel("$am_{\pi}$")
plt.show()