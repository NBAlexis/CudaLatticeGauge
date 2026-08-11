import lmfit
from lmfit import Model, Parameter, report_fit
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, fsolve

from CorrelatedFit import correlated_fit_with_flexible_bounds
from MesonStructures import all_signs, all_deltas, deltaidx_to_delta, delta_to_deltaidx

nt = 48
conf_count = 51
all_signs_lst = all_signs()
all_delta_lst = all_deltas()


def my_func(t, args):
    Ap = args[0]
    mp = args[1]
    Am = args[2]
    mm = args[3]
    return Ap * (np.exp(-mp * t) + np.exp(-mp * (48 - t))) + ((-1)**t) * Am * (np.exp(-mm * t) + np.exp(-mm * (48 - t)))

correlation_func = np.zeros(nt, dtype=np.complex128)
for type in range(20):
    all_type = len(all_signs_lst[type])
    for j in range(all_type):
        sign_func = all_signs_lst[type][j]
        delta_func = all_delta_lst[type][j]
        d = delta_to_deltaidx(delta_func)
        sum_res_lst = np.zeros((conf_count, nt), dtype=np.complex128)
        print(f"===={type}-{j}====")
        for conf in range(conf_count):
            p_nt_a_b_d = np.load(f"data/p_{1000 + 20 * conf}.npy")
            for t in range(nt):
                # find s(A) s(B) and delta to do the sum
                for av_idx in range(8):
                    for bv_idx in range(8):
                        av = deltaidx_to_delta(av_idx)
                        bv = deltaidx_to_delta(bv_idx)
                        sa = np.array([sign_func[m] * av[m] for m in range(3)])
                        sb = np.array([sign_func[m] * bv[m] for m in range(3)])
                        sign_prod = (-1) ** (np.sum(sa) + np.sum(sb))
                        sum_res_lst[conf, t] = sum_res_lst[conf, t] + sign_prod * p_nt_a_b_d[t, av_idx, bv_idx, d]
        print(f"saving : {type} - {j} : sn = {all_signs_lst[type][j]} delta={all_delta_lst[type][j]}")
        np.save(f"data/correlationpp_{type}_{j}.npy", sum_res_lst)
