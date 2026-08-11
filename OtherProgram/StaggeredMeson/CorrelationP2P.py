"""
for each sign_function, and delta

C(t) = sum_{x,y1,y2,A,B,c1,c2} s(A)s(B) D^-1_{2x+A,2y1+B;c1,c2} Dd^-1*_{2x+A+d,2y2+B+d;c1,c2}

since: s1_{B,c2}(x,c1) = sum_{y} D^-1(x, 2y+B; c1, c2)
       s2_{B,c2}(x,c1) = sum_{y} Dd^-1(x, 2y+B; c1, c2)

       so
       sum_{x,y1,y2,c1,c2} D^-1_{2x+A,2y1+B;c1,c2} Dd^-1*_{2x+A+d,2y2+B+d;c1,c2}
        = sum_{x,c1} s1_{B,c2}(2x+A,c1) s2*_{B+d,c2}(2x+A+d,c1)

        so, we pick 2x+A from s1_{B,c2} and 2x+A+d from s2_{B+d,c2}

        and then do the dot (sum over x)


    The quantity to be saved is p_{A,B,d}(t) = sum_{x,y1,y2,c1,c2} D^-1_{2x+A,2y1+B;c1,c2} Dd^-1*_{2x+A+d,2y2+B+d;c1,c2}
    and C(t) = sum _{A,B} s(A) s(B) p_{A,B,d}(t)

"""
import numpy as np

from MesonStructures import deltaidx_to_delta, delta_to_deltaidx, all_signs, all_deltas, third_sign

# 可选：若安装了 numba 可取消下面两行的注释以获得额外加速
# from numba import njit
# @njit
def dot_over_space_color_batch(pick1, pick2):
    """
    pick1, pick2 形状: (nt, A, nx//2, ny//2, nz//2, 3)
    返回形状: (nt, A)
    """
    return np.einsum('taxyzc,taxyzc->ta', pick1, pick2)

foldername1 = "G:\\mass"
foldername2 = "e32p31"
filenames = "propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371"

m = 0.00743 * 2
nx = 32
nt = 64
conf_count = 51

p_nt_a_b_d = np.zeros((conf_count, nt, 8, 8, 8), dtype=np.complex128)

# """
for conf in range(conf_count):
    conf_idx = 1000 + 20 * conf
    print("calculating conf = " + str(conf_idx))
    for A in range(8):
        for B in range(8):
            for delta in range(8):
                p_nt_a_b_d[conf, :, A, B, delta] = 0
                for c1 in range(3):
                    # propagator1 is pick from B+d and sink at A+d
                    # propagator2 is pick from B and sink at A
                    av = deltaidx_to_delta(A)
                    bv = deltaidx_to_delta(B)
                    dv = deltaidx_to_delta(delta)
                    avd = (dv + av) % 2
                    bvd = (dv + bv) % 2
                    avd_idx = delta_to_deltaidx(avd)
                    bvd_idx = delta_to_deltaidx(bvd)
                    psi1 = np.load(
                        f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{bvd_idx}_c{c1 + 1}_lexico.npy")
                    psi2 = np.load(
                        f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{bvd_idx}_dc{c1 + 1}_lexico.npy")
                    phi1 = -psi2 + m * psi1
                    psi3 = np.load(
                        f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{B}_c{c1 + 1}_lexico.npy")
                    psi4 = np.load(
                        f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{B}_dc{c1 + 1}_lexico.npy")
                    phi2 = psi4 + m * psi3
                    # now pick a+d from phi1 and a from phi2
                    pick1 = phi1[:, ((avd_idx >> 2) & 1):nx:2, ((avd_idx >> 1) & 1):nx:2, (avd_idx & 1):nx:2, :]
                    pick2 = phi2[:, ((A >> 2) & 1):nx:2, ((A >> 1) & 1):nx:2, (A & 1):nx:2, :]
                    # the shape of pick1 and pick2 should be the same, so just multiply them
                    pick1 = np.conj(pick1)
                    res = pick1 * pick2
                    # sum over x,y,z and c2
                    res = np.sum(res, axis=(1, 2, 3, 4))
                    # res is res[t]
                    p_nt_a_b_d[conf, :, A, B, delta] = p_nt_a_b_d[conf, :, A, B, delta] + res

# """


"""

# ---------- 第一组大循环（优化版） ----------
for conf in range(conf_count):
    conf_idx = 1000 + 20 * conf
    print("calculating conf = " + str(conf_idx))

    for c1 in range(3):
        # ---------- 缓存 1：仅依赖 (B, c1) 的数据 ----------
        phi2_parity_cache = {}   # key = B, value = phi2 的 8 种奇偶切片 (nt, 8, nx//2, ny//2, nz//2, 3)

        for B in range(8):
            psi3 = np.load(f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{B}_c{c1 + 1}_lexico.npy")
            psi4 = np.load(f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{B}_dc{c1 + 1}_lexico.npy")
            phi2 = psi4 + m * psi3

            # 预先提取 phi2 的 8 种奇偶起始点的子格点（对应 8 个 A）
            slices = []
            for A in range(8):
                x_start = (A >> 2) & 1
                y_start = (A >> 1) & 1
                z_start = A & 1
                slices.append(phi2[:, x_start::2, y_start::2, z_start::2, :])
            phi2_parity_cache[B] = np.stack(slices, axis=1)   # (nt, 8, nx//2, ny//2, nz//2, 3)

        # ---------- 缓存 2：仅依赖 (B, delta, c1) 的数据 ----------
        phi1_parity_cache = {}   # key = (B, delta), value = phi1 的 8 种奇偶切片（共轭后）

        for B in range(8):
            bv = deltaidx_to_delta(B)
            for delta in range(8):
                dv = deltaidx_to_delta(delta)
                bvd = (dv + bv) % 2
                bvd_idx = delta_to_deltaidx(bvd)

                psi1 = np.load(f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{bvd_idx}_c{c1 + 1}_lexico.npy")
                psi2 = np.load(f"{foldername1}\\{foldername2}\\{filenames}_{conf_idx}_delta{bvd_idx}_dc{c1 + 1}_lexico.npy")
                phi1 = -psi2 + m * psi1

                slices = []
                for A in range(8):
                    # 注意：phi1 选取的是 A^delta 对应的奇偶偏移
                    avd_idx = A ^ delta   # 等价于 (A + delta) % 8 的位运算形式
                    x_start = (avd_idx >> 2) & 1
                    y_start = (avd_idx >> 1) & 1
                    z_start = avd_idx & 1
                    slices.append(phi1[:, x_start::2, y_start::2, z_start::2, :])
                phi1_parity_cache[(B, delta)] = np.conj(np.stack(slices, axis=1))

        # ---------- 计算并累加（此时消除了 A 循环） ----------
        for B in range(8):
            for delta in range(8):
                pick2 = phi2_parity_cache[B]                 # (nt, 8, nx//2, ny//2, nz//2, 3)
                pick1 = phi1_parity_cache[(B, delta)]        # (nt, 8, nx//2, ny//2, nz//2, 3)
                res_nt_A = dot_over_space_color_batch(pick1, pick2)   # (nt, 8)
                p_nt_a_b_d[conf, :, :, B, delta] += res_nt_A


"""

all_signs_lst = all_signs()
all_delta_lst = all_deltas()

correlation_func = np.zeros(nt, dtype=np.complex128)
for type in range(20):
    all_type = len(all_signs_lst[type])
    for j in range(all_type):
        sign_func = all_signs_lst[type][j]
        delta_func = all_delta_lst[type][j]
        # some k=3, l=1, m=2 channels have different sign
        # third_channel_sign = 1 if j != 2 else all_third_channel_phase[type]
        d = delta_to_deltaidx(delta_func)
        sum_res_lst = np.zeros((conf_count, nt), dtype=np.complex128)
        print(f"===={type}-{j}====")
        for conf in range(conf_count):
            p_nt_a_b_d_this_conf = p_nt_a_b_d[conf, :, :, :, :]
            for t in range(nt):
                # find s(A) s(B) and delta to do the sum
                for av_idx in range(8):
                    for bv_idx in range(8):
                        av = deltaidx_to_delta(av_idx)
                        bv = deltaidx_to_delta(bv_idx)
                        sa = np.array([sign_func[m] * av[m] for m in range(3)])
                        sb = np.array([sign_func[m] * bv[m] for m in range(3)])
                        sign_prod = (-1) ** (np.sum(sa) + np.sum(sb))
                        sum_res_lst[conf, t] = sum_res_lst[conf, t] + sign_prod * p_nt_a_b_d_this_conf[t, av_idx, bv_idx, d]
        print(f"saving : {type} - {j} : sn = {all_signs_lst[type][j]} delta={all_delta_lst[type][j]}")
        np.save(f"data/{foldername2}/correlationp2p_{type}_{j}.npy", sum_res_lst)
