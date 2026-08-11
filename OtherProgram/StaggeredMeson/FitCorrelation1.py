import lmfit
from lmfit import Model, Parameter, report_fit
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, fsolve

from CorrelatedFit import correlated_fit_with_flexible_bounds
from CorrelatedFitLMFIT import  correlated_fit_with_flexible_bounds_lmfit
from MesonStructures import all_signs

nt = 48
t_start = 6
t_end = 23

for type in range(20):
    # type = 4
    all_signs_lst = all_signs()
    all_type = len(all_signs_lst[type])

    def my_func(t, args):
        Ap = args[0]
        mp = args[1]
        Am = args[2]
        mm = args[3]
        return Ap * (np.exp(-mp * t) + np.exp(-mp * (48 - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (48 - t)))

    correlation_func = np.load(f"data/correlation_{type}_0.npy")
    for i in range(1, all_type):
        correlation_func = correlation_func + np.load(f"data/correlation_{type}_{i}.npy")

    correlation_func = np.real(correlation_func)
    # pick a separation of t
    t_value = np.array([k for k in range(nt)])
    # t_stride = t_end - t_start
    # pickidx = [k + t_start + 1 for k in range(t_stride)] + [nt - t_stride - t_start + k for k in range(t_stride)]
    pickidx = [k + t_start + 1 for k in range(nt - 1 - 2 * t_start)]
    x_pick = t_value[pickidx]
    y_pick = correlation_func[:, pickidx]
    initial_guess = np.array([np.mean(correlation_func), 0.2, np.std(correlation_func), 0.2])
    upper_bounds = [np.inf, 1.0, np.inf, 1.0]
    lower_bounds = [-np.inf, 0.1, -np.inf, 0.1]

    """
        method : str
            lmfit 优化方法，应该与scipy的curve_fit对应：
            - 'least_squares': 最小二乘法（默认，对应'trf'/'dogbox'）
            - 'leastsq': Levenberg-Marquardt（对应'lm'）
            - 'lbfgsb': L-BFGS-B 算法
            - 'differential_evolution': 差分进化法（全局优化）
            - 'nelder': Nelder-Mead 单纯形法
    """

    fitres = correlated_fit_with_flexible_bounds_lmfit(
            y_pick,
            x_pick,
            my_func,
            initial_guess,
            param_bounds=(lower_bounds, upper_bounds)
    )
    print(f"channel:{type}, m+ {fitres['params'][1]} +- {fitres['params_err'][1]} ({fitres['params'][1]*1827.1} +- {fitres['params_err'][1]*1827.1}) and m- = {fitres['params'][3]} +- {fitres['params_err'][3]} ({fitres['params'][3]*1827.1} +- {fitres['params_err'][3]*1827.1}), chi2/dof = {fitres['chi2_dof']}")

