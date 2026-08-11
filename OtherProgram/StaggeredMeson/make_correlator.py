"""
从已有的 p_{conf}.npy 文件生成 correlator 文件。
完全沿用 propav2.py 后半部分的逻辑，只是把第一步"QUDA 求逆 + 构建 p 数组"
替换为"读取 p_{conf}.npy"。

Usage: python make_correlator.py <folder_name>
  e.g. python make_correlator.py c48p13
  从 data/newcorrelator/c48p13/p_*.npy 读取，输出到同目录。
"""
import numpy as np
import sys
import os
import glob


# ── helpers from propav2.py ──

def eta_i(i):
    return np.array([1 if (j + 1) < i else 0 for j in range(4)])

def xi_i(i):
    return np.array([1 if (j + 1) > i else 0 for j in range(4)])

def spatial_sign(arr):
    return np.array([arr[0], arr[1], arr[2]])

def deltaidx_to_delta(delta):
    return np.array([delta & 1, (delta >> 1) & 1, (delta >> 2) & 1])

def delta_to_deltaidx(delta):
    return (delta[2] << 2) | (delta[1] << 1) | delta[0]

def all_signs():
    S = []
    S.append([np.array([0, 0, 0])])
    S.append([spatial_sign((eta_i(4) + xi_i(4)) % 2)])
    S.append([spatial_sign((eta_i(k) + xi_i(k) + eta_i(5)) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(k) + xi_i(k) + eta_i(5)) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign(eta_i(k) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(k)) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign((xi_i(k) + eta_i(5)) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + xi_i(k) + eta_i(5)) % 2) for k in (1, 2, 3)])
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(j) + eta_i(5)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(j) + eta_i(5)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(i) + eta_i(j)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + eta_i(j)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((xi_i(i) + xi_i(j)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + xi_i(i) + xi_i(j)) % 2)
              for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(j) + xi_i(k)) % 2)
              for i, j, k in [(3, 1, 2), (1, 2, 3), (2, 3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(j) + xi_i(k)) % 2)
              for i, j, k in [(3, 1, 2), (1, 2, 3), (2, 3, 1)]])
    S.append([spatial_sign((eta_i(1) + eta_i(2) + eta_i(3)) % 2)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + eta_i(2) + eta_i(3)) % 2)])
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3)) % 2)
              for i in (1, 2, 3)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3)) % 2)
              for i in (1, 2, 3)])
    return S

def all_deltas():
    D = []
    D.append([np.array([0, 0, 0])])
    D.append([np.array([0, 0, 0])])
    D.append([np.array([0, 0, 0])] * 3)
    D.append([np.array([0, 0, 0])] * 3)
    D.append([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    D.append([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    D.append([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    D.append([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    D.append([np.array([0, 1, 0]), np.array([0, 0, 1]), np.array([1, 0, 0])])
    D.append([np.array([0, 1, 0]), np.array([0, 0, 1]), np.array([1, 0, 0])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 0]), np.array([0, 1, 1]), np.array([1, 0, 1])])
    D.append([np.array([1, 1, 1])])
    D.append([np.array([1, 1, 1])])
    D.append([np.array([1, 1, 1])] * 3)
    D.append([np.array([1, 1, 1])] * 3)
    return D


def generate_correlators(in_dir, out_dir):
    """Read p_*.npy from in_dir, write correlator files to out_dir."""
    p_files = sorted(glob.glob(f"{in_dir}/p_*.npy"),
                     key=lambda x: int(x.split('_')[-1].replace('.npy', '')))
    conf_count = len(p_files)
    if conf_count == 0:
        raise FileNotFoundError(f"No p_*.npy files found in {in_dir}")

    first = np.load(p_files[0])
    nt = first.shape[0]
    print(f"[{in_dir}] Loading {conf_count} p files, nt={nt}")

    p_nt_a_b_d = np.zeros((conf_count, nt, 8, 8, 8), dtype=np.complex128)
    for conf, path in enumerate(p_files):
        p_nt_a_b_d[conf] = np.load(path)

    all_signs_lst = all_signs()
    all_delta_lst = all_deltas()
    basestr = os.path.basename(in_dir)

    os.makedirs(out_dir, exist_ok=True)

    for type in range(20):
        all_type = len(all_signs_lst[type])
        for j in range(all_type):
            sign_func = all_signs_lst[type][j]
            delta_func = all_delta_lst[type][j]
            d = delta_to_deltaidx(delta_func)
            sum_res_lst = np.zeros((conf_count, nt), dtype=np.complex128)
            for conf in range(conf_count):
                p_conf = p_nt_a_b_d[conf]
                for t in range(nt):
                    for av_idx in range(8):
                        for bv_idx in range(8):
                            av = deltaidx_to_delta(av_idx)
                            bv = deltaidx_to_delta(bv_idx)
                            sa = np.array([sign_func[m] * av[m] for m in range(3)])
                            sb = np.array([sign_func[m] * bv[m] for m in range(3)])
                            sign_prod = (-1) ** (np.sum(sa) + np.sum(sb))
                            sum_res_lst[conf, t] += sign_prod * p_conf[t, av_idx, bv_idx, d]
            np.save(f"{out_dir}/correlationp2p_{type}_{j}.npy", sum_res_lst)

    for channel in range(20):
        combined = np.load(f"{out_dir}/correlationp2p_{channel}_0.npy")
        for i in range(1, len(all_signs_lst[channel])):
            combined += np.load(f"{out_dir}/correlationp2p_{channel}_{i}.npy")
        np.save(f"{out_dir}/{basestr}_{channel}.npy", np.real(combined))


# ── main ──
if __name__ == "__main__":
    import sys

    base = "data/newcorrelator"
    if len(sys.argv) > 1:
        # 指定文件夹: python make_correlator.py x24p31 z24p31
        targets = sys.argv[1:]
    else:
        # 无参数: 处理所有含 p_*.npy 的文件夹
        targets = sorted(
            os.path.basename(d.rstrip("/\\"))
            for d in glob.glob(f"{base}/*/")
            if glob.glob(os.path.join(d, "p_*.npy"))
        )

    print(f"Folders to process: {targets}\n")
    for basestr in targets:
        d = f"{base}/{basestr}"
        if not glob.glob(f"{d}/p_*.npy"):
            print(f"=== SKIP {basestr}: no p_*.npy ===")
            continue
        print(f"=== Processing {basestr} ===")
        generate_correlators(d, d)
        p2p_cnt = len(glob.glob(f"{d}/correlationp2p_*.npy"))
        comb_cnt = len(glob.glob(f"{d}/{basestr}_*.npy"))
        print(f"  Done: {p2p_cnt} correlationp2p + {comb_cnt} combined files\n")
