"""
propav3_mpi.py  (基于 propav3.py，适配 4 进程 MPI + 4 DCU)
=============================================================
每个 MPI 进程控制 1 块 DCU，QUDA 内部通过 MPI 通信并行求解。
I/O 和最终分析只由 rank 0 执行，避免文件冲突。
"""
import numpy as np
import sys
import os

# ---------- MPI 初始化 ----------
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# import os
# os.environ['ROCR_VISIBLE_DEVICES'] = str(rank % 4)
# os.environ['HIP_VISIBLE_DEVICES'] = str(rank % 4)
# os.environ['CUDA_VISIBLE_DEVICES'] = str(rank % 4)

# ---------- QUDA imports ----------
from pyquda import init
from pyquda_utils import core, io
from pyquda_utils.io import readNPYGauge, writeNPYGauge
from pyquda_utils.core import LatticeInfo
from pyquda.field import LatticeGauge
from pyquda.field import (
    LatticeStaggeredFermion,
    MultiLatticeStaggeredFermion,
    LatticeStaggeredPropagator,
)
from pyquda.enum_quda import QudaSolutionType, QudaSolveType, QudaInverterType

# ---------- 参数解析 ----------
args = sys.argv[1:]
nx = int(args[0])
nt = int(args[1])
mass = float(args[2])
folder = args[3]
filename = args[4]
basestr = args[5]
nums = [int(args[6]), int(args[7]), int(args[8])]

if rank == 0:
    print(f"Running with {size} MPI processes")
    print(f"nx={nx}  nt={nt}  mass={mass}")
    print(f"folder={folder}")
    print(f"filename={filename}")
    print(f"basestr={basestr}")
    print(f"start_idx = {nums[0]}  count = {nums[1]}  stride = {nums[2]}")

# ---------- Helper functions (unchanged) ----------
def ColorVectorWithDelta(latt_info, t, Nx, Ny, Nz, delta, color):
    Lx, Ly, Lz, Lt = latt_info.Lx, latt_info.Ly, latt_info.Lz, latt_info.Lt
    gx, gy, gz, gt = latt_info.gx, latt_info.gy, latt_info.gz, latt_info.gt
    dx, dy, dz = delta & 1, (delta >> 1) & 1, (delta >> 2) & 1
    b = LatticeStaggeredFermion(latt_info)
    if gt * Lt <= t < (gt + 1) * Lt:
        for x in range(dx, Nx, 2):
            for y in range(dy, Ny, 2):
                for z in range(dz, Nz, 2):
                    if gx * Lx <= x < (gx + 1) * Lx and gy * Ly <= y < (gy + 1) * Ly and gz * Lz <= z < (gz + 1) * Lz:
                        eo = ((x - gx * Lx) + (y - gy * Ly) + (z - gz * Lz) + (t - gt * Lt)) % 2
                        b.data[eo, t - gt * Lt, z - gz * Lz, y - gy * Ly, (x - gx * Lx) // 2, color] = 1
    return b

def invertStaggeredWithDelta(dirac, t, delta):
    latt_info = dirac.latt_info
    propag = LatticeStaggeredPropagator(latt_info)
    for color in range(3):
        b = MultiLatticeStaggeredFermion(latt_info, 1)
        b[0] = ColorVectorWithDelta(latt_info, t, nx, nx, nx, delta, color)
        x = dirac.invertMultiSrcRestart(b, 0)
        propag.setFermion(x[0], color)
    return propag

# ---------- Meson structure functions (unchanged) ----------
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
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(j) + eta_i(5)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(j) + eta_i(5)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(i) + eta_i(j)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + eta_i(j)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((xi_i(i) + xi_i(j)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + xi_i(i) + xi_i(j)) % 2) for i, j in [(1, 2), (2, 3), (3, 1)]])
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(j) + xi_i(k)) % 2) for i, j, k in [(3, 1, 2), (1, 2, 3), (2, 3, 1)]])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(j) + xi_i(k)) % 2) for i, j, k in [(3, 1, 2), (1, 2, 3), (2, 3, 1)]])
    S.append([spatial_sign((eta_i(1) + eta_i(2) + eta_i(3)) % 2)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + eta_i(2) + eta_i(3)) % 2)])
    S.append([spatial_sign((eta_i(i) + xi_i(i) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3)) % 2) for i in (1, 2, 3)])
    S.append([spatial_sign((eta_i(4) + xi_i(4) + eta_i(i) + xi_i(i) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3)) % 2) for i in (1, 2, 3)])
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

# ---------- QUDA 初始化 ----------
# grid_size 乘积必须等于 MPI 进程数 (这里是 1*1*2*2 = 4)
init(grid_size=[2, 2, 2, 1], latt_size=[nx, nx, nx, nt], resource_path=".cache", backend="numpy")
latt_info = core.LatticeInfo([nx, nx, nx, nt], -1, 1.0)

out_dir = f"data/{basestr}"
if rank == 0:
    os.makedirs(out_dir, exist_ok=True)

m = mass * 2
conf_count = nums[1]

all_signs_lst = all_signs()
all_delta_lst = all_deltas()

# 只在 rank 0 分配完整的结果数组，其他进程 None
if rank == 0:
    p_nt_a_b_d = np.zeros((conf_count, nt, 8, 8, 8), dtype=np.complex128)
else:
    p_nt_a_b_d = None

for conf in range(conf_count):
    confidx = nums[0] + nums[2] * conf

    # 所有进程一同读取 gauge（避免用 broadcast 复杂化）
    gauge = io.readChromaQIOGauge(f'{folder}/{filename}/{filename}_{confidx}_coulomb_1e-9.scidac')
    if rank == 0:
        print(f"gauge loaded... conf={confidx}")

    dslash = core.getHISQ(latt_info, mass, 1e-12, 10000, 0.0)
    dslash.invert_param.solution_type = QudaSolutionType.QUDA_MATDAG_MAT_SOLUTION
    dslash.invert_param.solve_type = QudaSolveType.QUDA_NORMOP_SOLVE
    dslash0 = core.getHISQ(latt_info, 0.0, 1e-12, 10000, 0.0)        

    dslash.loadGauge(gauge)
    dslash0.loadGauge(gauge)

    fermion_delta_c1, fermion_delta_c2, fermion_delta_c3 = [], [], []
    dfermion_delta_c1, dfermion_delta_c2, dfermion_delta_c3 = [], [], []

    for delta in range(8):
        propagator = invertStaggeredWithDelta(dslash, 0, delta)
        fermion_0 = propagator.getFermion(0)
        fermion_1 = propagator.getFermion(1)
        fermion_2 = propagator.getFermion(2)

        fermion_delta_c1.append(np.asarray(fermion_0.lexico()))
        fermion_delta_c2.append(np.asarray(fermion_1.lexico()))
        fermion_delta_c3.append(np.asarray(fermion_2.lexico()))

        dfermion_0 = dslash0.dslash(fermion_0)
        dfermion_1 = dslash0.dslash(fermion_1)
        dfermion_2 = dslash0.dslash(fermion_2)
        dfermion_delta_c1.append(np.asarray(dfermion_0.lexico()))
        dfermion_delta_c2.append(np.asarray(dfermion_1.lexico()))
        dfermion_delta_c3.append(np.asarray(dfermion_2.lexico()))

    fermion_delta_c = [fermion_delta_c1, fermion_delta_c2, fermion_delta_c3]
    dfermion_delta_c = [dfermion_delta_c1, dfermion_delta_c2, dfermion_delta_c3]

    # 每个进程只在本地格点部分计算贡献
    local_p_slice = np.zeros((nt, 8, 8, 8), dtype=np.complex128)

    for A in range(8):
        for B in range(8):
            for delta in range(8):
                local_sum_t = np.zeros(nt, dtype=np.complex128)
                for c1 in range(3):
                    av = deltaidx_to_delta(A)
                    bv = deltaidx_to_delta(B)
                    dv = deltaidx_to_delta(delta)
                    avd = (dv + av) % 2
                    bvd = (dv + bv) % 2
                    avd_idx = delta_to_deltaidx(avd)
                    bvd_idx = delta_to_deltaidx(bvd)

                    psi1 = fermion_delta_c[c1][bvd_idx]
                    psi2 = dfermion_delta_c[c1][bvd_idx]
                    phi1 = -psi2 + m * psi1

                    psi3 = fermion_delta_c[c1][B]
                    psi4 = dfermion_delta_c[c1][B]
                    phi2 = psi3 + m * psi4

                    pick1 = phi1[:, ((avd_idx >> 2) & 1):nx:2,
                                    ((avd_idx >> 1) & 1):nx:2,
                                    (avd_idx & 1):nx:2, :]
                    pick2 = phi2[:, ((A >> 2) & 1):nx:2,
                                    ((A >> 1) & 1):nx:2,
                                    (A & 1):nx:2, :]
                    res = np.conj(pick1) * pick2
                    local_sum_t += np.sum(res, axis=(1, 2, 3, 4))

                # 把各进程的局部贡献加起来（Allreduce），得到全局的 nt 向量
                comm.Allreduce(MPI.IN_PLACE, local_sum_t, op=MPI.SUM)
                local_p_slice[:, A, B, delta] = local_sum_t

    # rank 0 收集结果并保存
    if rank == 0:
        p_nt_a_b_d[conf] = local_p_slice
        np.save(f"{out_dir}/p_{conf}.npy", p_nt_a_b_d[conf])
        print(f"  correlators conf={confidx} done")

    del dslash, dslash0, fermion_delta_c, dfermion_delta_c
    # cupy.cuda.Device().synchronize()
    # mempool = cupy.get_default_memory_pool()
    # mempool.free_all_blocks()
    # pinned_mempool = cupy.get_default_pinned_memory_pool()
    # pinned_mempool.free_all_blocks()    

# ---------- 投影到介子通道 (仅 rank 0) ----------
if rank == 0:
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
                p_this = p_nt_a_b_d[conf]
                for t in range(nt):
                    for av_idx in range(8):
                        for bv_idx in range(8):
                            av = deltaidx_to_delta(av_idx)
                            bv = deltaidx_to_delta(bv_idx)
                            sa = np.array([sign_func[m] * av[m] for m in range(3)])
                            sb = np.array([sign_func[m] * bv[m] for m in range(3)])
                            sign_prod = (-1) ** (np.sum(sa) + np.sum(sb))
                            sum_res_lst[conf, t] += sign_prod * p_this[t, av_idx, bv_idx, d]
            print(f"saving : {type} - {j} : sn = {all_signs_lst[type][j]} delta={all_delta_lst[type][j]}")
            np.save(f"{out_dir}/correlationp2p_{type}_{j}.npy", sum_res_lst)

    for channel in range(20):
        correlation_func_combine = np.load(f"{out_dir}/correlationp2p_{channel}_0.npy")
        for i in range(1, len(all_signs_lst[channel])):
            correlation_func_combine += np.load(f"{out_dir}/correlationp2p_{channel}_{i}.npy")
        np.save(f"{out_dir}/{basestr}_{channel}.npy", np.real(correlation_func_combine))

    print("All done.")