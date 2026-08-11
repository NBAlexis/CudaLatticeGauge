"""
for each sign_function, and delta

C(t) = sum_{A,B,c1,c2} s(A)s(B) u_{A,B;c1,c2} v*_{A+d,B+d;c1,c2}

"""
import numpy as np

from MesonStructures import all_signs, all_deltas, delta_to_deltaidx, deltaidx_to_delta

all_signs_lst = all_signs()
all_delta_lst = all_deltas()
nt = 48

ct_all = []

beforesum = np.zeros((nt, 8, 8, 3, 3), dtype=np.complex128)
for conf in range(51):
    conf_idx = 1000 + 20 * conf
    ct_this_configuration = []
    print("calculating conf = " + str(conf_idx))
    u_times_v = np.load(f"data/u{conf_idx}.npy")
    v = np.load(f"data/v{conf_idx}.npy")
    for i in range(20):
        all_type = len(all_signs_lst[i])
        for j in range(all_type):
            print(f"===={i}-{j}====")
            sign_func = all_signs_lst[i][j]
            delta_func = all_delta_lst[i][j]
            sum_res_lst = np.zeros(nt, dtype=np.complex128)
            for t in range(nt):
                sum_res = 0j
                for av_idx in range(8):
                    for bv_idx in range(8):
                        # calculate av+d and bv+d
                        av = deltaidx_to_delta(av_idx)
                        bv = deltaidx_to_delta(bv_idx)
                        avd = (delta_func + av)%2
                        bvd = (delta_func + bv)%2
                        avd_idx = delta_to_deltaidx(avd)
                        bvd_idx = delta_to_deltaidx(bvd)
                        sa = np.array([sign_func[m] * av[m] for m in range(3)])
                        sb = np.array([sign_func[m] * bv[m] for m in range(3)])
                        sign_prod = (-1) ** (np.sum(sa) + np.sum(sb))
                        for c1 in range(3):
                            for c2 in range(3):
                                """
                                u[t, a, b, c1, c2] = sink[a,c1].D-1.source[b,c2]
                                """
                                # sum_res = sum_res + sign_prod * u[t, av_idx, bv_idx, c1, c2] * np.conj(v[t, bvd_idx, avd_idx, c2, c1])
                                sum_res = sum_res + sign_prod * u[t, av_idx, bv_idx, c1, c2] * np.conj(v[t, avd_idx, bvd_idx, c1, c2])
                sum_res_lst[t] = sum_res
            ct_this_configuration.append(sum_res_lst)
    ct_all.append(ct_this_configuration)
ct_all = np.array(ct_all)
print(np.shape(ct_all))

k = 0
for i in range(20):
    all_type = len(all_signs_lst[i])
    for j in range(all_type):
        # ct_all[configuration, type, t]
        print(f"saving : {i} - {j} : sn = {all_signs_lst[i][j]} delta={all_delta_lst[i][j]}")
        np.save(f"data/correlation_{i}_{j}.npy", ct_all[:,k,:])
        k = k + 1




