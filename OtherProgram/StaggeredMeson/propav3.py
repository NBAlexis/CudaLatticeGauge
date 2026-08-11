"""
propav2.py
==========
Fused propagator-inversion + correlator-construction in a single pass.

For each configuration:
  1. Compute  w1 = (D+D)^{-1} w_delta   for all delta=0..7, color=0..2
  2. Compute  w2 = D0. 1
  3. Build phi1 = -w2 + m*w1,  phi2 = w1 + m*w2
  4. Contract into p_{A,B,delta}(t) and project onto all 20 meson channels

Propagators are kept in memory only - no intermediate files written.
Output: one .npy file per (type, j) channel, same format as CorrelationP2P.py.


cd /public/home/nbalexis/pyquda/working
source ./env.sh
python propav2.py 24 48 0.00944 /public/share/ybyang/caimc/CloveronHISQ s24t48_beta7.29_ml0.00944ms0.04721mc0.5555 c24p31 1000 51 20
python propav2.py 24 48 0.00944 /public/share/ybyang/ltw/CloveronHISQ s24t48_beta7.29_ml0.00944ms0.055mc0.5555 c24p31s 1000 51 20
python propav2.py 24 48 0.00472 /public/share/ybyang/ltw/CloveronHISQ s24t48_beta7.29_ml0.00472ms0.04721mc0.5555 c24p22 1000 51 20
python propav2.py 32 48 0.00944 /public/share/ybyang/ltw/CloveronHISQ s32t48_beta7.29_ml0.00944ms0.04721mc0.5555 c32p31 1000 51 20
python propav2.py 32 48 0.00472 /public/share/ybyang/ltw/CloveronHISQ s32t48_beta7.29_ml0.00472ms0.04721mc0.5555 c32p22 1000 51 20
python propav2.py 48 48 0.00174 /public/share/ybyang/ltw/CloveronHISQ s48t48_beta7.29_ml0.00174ms0.04721mc0.5555 c48p13 1000 51 20
python propav2.py 32 64 0.00743 /public/share/ybyang/caimc/CloveronHISQ s32t64_beta7.54_ml0.00743ms0.03715mc0.4371 e32p31 1000 51 20
python propav2.py 32 64 0.00579 /public/share/ybyang/caimc/CloveronHISQ s32t64_beta7.75_ml0.00579ms0.02895mc0.34 g32p32 1000 51 20
python propav2.py 48 64 0.00579 /public/share/ybyang/ltw/CloveronHISQ s48t64_beta7.75_ml0.00579ms0.02895mc0.34 g48p31 1000 51 20
python propav2.py 48 96 0.00352 /public/share/ybyang/caimc/CloveronHISQ s48t96_beta8.2_ml0.003526ms0.01763mc0.207 h48p31 1000 51 20
python propav2.py 24 48 0.009796 /public/share/ybyang/ltw/CloveronHISQ s24t48_beta7.2133_ml0.009796ms0.04898 y24p31 1000 51 20



c24P31
/public/share/ybyang/caimc/CloveronHISQ/s24t48_beta7.29_ml0.00944ms0.04721mc0.5555/s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf}_coulomb_1e-9.scidac
c24P31s
/public/share/ybyang/ltw/CloveronHISQ/s24t48_beta7.29_ml0.00944ms0.055mc0.5555/s24t48_beta7.29_ml0.00944ms0.055mc0.5555_{conf}_coulomb_1e-9.scidac
c24P22
/public/share/ybyang/ltw/CloveronHISQ/s24t48_beta7.29_ml0.00472ms0.04721mc0.5555/s24t48_beta7.29_ml0.00472ms0.04721mc0.5555_{conf}_coulomb_1e-9.scidac
c32P31
/public/share/ybyang/ltw/CloveronHISQ/s32t48_beta7.29_ml0.00944ms0.04721mc0.5555/s32t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf}_coulomb_1e-9.scidac
c32P22
/public/share/ybyang/ltw/CloveronHISQ/s32t48_beta7.29_ml0.00472ms0.04721mc0.5555/s32t48_beta7.29_ml0.00472ms0.04721mc0.5555_{conf}_coulomb_1e-9.scidac
c48P13
/public/share/ybyang/ltw/CloveronHISQ/s48t48_beta7.29_ml0.00174ms0.04721mc0.5555/s48t48_beta7.29_ml0.00174ms0.04721mc0.5555_{conf}_coulomb_1e-9.scidac
e32P31
/public/share/ybyang/caimc/CloveronHISQ/s32t64_beta7.54_ml0.00743ms0.03715mc0.4371/s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{conf}_coulomb_1e-9.scidac
g32P32
/public/share/ybyang/caimc/CloveronHISQ/s32t64_beta7.75_ml0.00579ms0.02895mc0.34/s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{conf}_coulomb_1e-9.scidac
g48P31
/public/share/ybyang/ltw/CloveronHISQ/s48t64_beta7.75_ml0.00579ms0.02895mc0.34/s48t64_beta7.75_ml0.00579ms0.02895mc0.34_{conf}_coulomb_1e-9.scidac
h48P31
/public/share/ybyang/caimc/CloveronHISQ/s48t96_beta8.2_ml0.003526ms0.01763mc0.207/s48t96_beta8.2_ml0.003526ms0.01763mc0.207_{conf}_coulomb_1e-9.scidac
y24P31
/public/share/ybyang/ltw/CloveronHISQ/s24t48_beta7.2133_ml0.009796ms0.04898/s24t48_beta7.2133_ml0.009796ms0.04898_{1000}_coulomb_1e-9.scidac

"""

# ── Lattice parameters from command line ──
# Usage: python propav2.py <nx> <nt> <mass> <folder> <filename> <basestr> <start_idx> <count> <stride>


# ── QUDA imports ──
import numpy as np
import sys
import os

from pyquda import init
from pyquda_utils import core, io
from pyquda_utils.io import readNPYGauge, writeNPYGauge
from pyquda_utils.core import LatticeInfo
from pyquda.field import (
    LatticeStaggeredFermion,
    MultiLatticeStaggeredFermion,
    LatticeStaggeredPropagator,
)
from pyquda.enum_quda import QudaSolutionType, QudaSolveType, QudaInverterType


args = sys.argv[1:]
nx = int(args[0])
nt = int(args[1])
mass = float(args[2])
folder = args[3]
filename = args[4]
basestr = args[5]
nums = [int(args[6]), int(args[7]), int(args[8])]
print(f"nx={nx}  nt={nt}  mass={mass}")
print(f"folder={folder}")
print(f"filename={filename}")
print(f"basestr={basestr}")
print(f"start_idx = {nums[0]}  count = {nums[1]}  stride = {nums[2]}")

# ──────────────────────────────────────────────────────────────────────────────
# Helper functions (ported from propa.py)
# ──────────────────────────────────────────────────────────────────────────────

def ColorVectorWithDelta(latt_info: LatticeInfo, t: int, Nx: int, Ny: int, Nz: int, delta: int, color: int):
    """
    Docstring for colorvector
    
    :param latt_info: Description
    :type latt_info: LatticeInfo
    :param t_srce: Description
    :type t_srce: int
    :param delta: 位移矢量,0到7，x,y,z轴是否为1分别为 delta&1, delta&2, delta&4
    :type delta: int
    :param color:
    :type color: int
    """
    Lx = latt_info.Lx
    Ly = latt_info.Ly
    Lz = latt_info.Lz
    Lt = latt_info.Lt
    gx = latt_info.gx
    gy = latt_info.gy
    gz = latt_info.gz
    gt = latt_info.gt
    # if 0 == delta and 0 == color:
    #     print(Lx)
    #     print(Ly)
    #     print(Lz)
    #     print(Lt)
    #     print(gx)
    #     print(gy)
    #     print(gz)
    #     print(gt)
    dx = delta&1
    dy = (delta>>1)&1
    dz = (delta>>2)&1
    """
    我们需要b.data[t, dz:Nz:2, dy:Ny:2. dx:Nx:2] = 1
    其中，even odd = (t + dx + dy + dz)&1
    考虑gx,有两种情况：
    当gx*Lx是偶数，那么从gx*Lx开始，否则从gx*Lx + 1开始
    """
    b = LatticeStaggeredFermion(latt_info)
    """
    这个可能是对的，不过下面的繁琐的代码虽然慢，也没有太慢，可以接受
    if gt * Lt <= t < (gt + 1) * Lt:
        even_odd = (t + dx + dy + dz) & 1
        x_start = gx * Lx
        y_start = gy * Ly
        z_start = gz * Lz
        x_start = x_start + ((x_start+gx)&1)
        y_start = y_start + ((y_start+gy)&1)
        z_start = z_start + ((z_start+gz)&1)
        x_end = (gx+1) * Lx
        y_end = (gy+1) * Ly
        z_end = (gz+1) * Lz        
        b.data[even_odd, t, z_start:z_end:2, y_start:y_end:2, x_start:x_end:2, color] = 1
    """
    if gt*Lt <= t < (gt+1)*Lt:
        for x in range(dx,Nx,2):
            for y in range(dy,Ny,2):
                for z in range(dz,Nz,2):
                    if gx*Lx <= x < (gx+1)*Lx:
                        if gy*Ly <= y < (gy+1)*Ly:
                            if gz*Lz <= z < (gz+1)*Lz:
                                eo = ((x - gx * Lx) + (y - gy * Ly) + (z - gz * Lz) + (t - gt * Lt)) % 2
                                b.data[eo, t - gt * Lt, z - gz * Lz, y - gy * Ly, (x - gx * Lx) // 2, color] = 1
    return b

def invertStaggeredWithDelta(dirac, t, delta):
    """
    我们就老老实实一个一个求解，L5=1
    
    :param dirac: Description
    :param t: Description
    :param delta: Description
    """
    latt_info = dirac.latt_info
    propag = LatticeStaggeredPropagator(latt_info)
    for color in range(3):
        b = MultiLatticeStaggeredFermion(latt_info, 1)
        b[0] = ColorVectorWithDelta(latt_info, t, nx, nx, nx, delta, color)
        # the second parameter is restart
        x = dirac.invertMultiSrcRestart(b, 0)
        propag.setFermion(x[0], color)
    return propag


# ──────────────────────────────────────────────────────────────────────────────
# Helper functions (ported from MesonStructures.py)
# ──────────────────────────────────────────────────────────────────────────────

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


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

init(grid_size=[1, 1, 2, 2], latt_size=[nx, nx, nx, nt], resource_path=".cache")
latt_info = core.LatticeInfo([nx, nx, nx, nt], -1, 1.0)

dslash = core.getHISQ(latt_info, mass, 1e-12, 5000, 0.0)
dslash.invert_param.solution_type = QudaSolutionType.QUDA_MATDAG_MAT_SOLUTION
dslash.invert_param.solve_type = QudaSolveType.QUDA_NORMOP_SOLVE

dslash0 = core.getHISQ(latt_info, 0.0, 1e-12, 5000, 0.0)

out_dir = f"data/{basestr}"
os.makedirs(out_dir, exist_ok=True)

m = mass * 2
conf_count = nums[1]

all_signs_lst = all_signs()
all_delta_lst = all_deltas()

p_nt_a_b_d = np.zeros((conf_count, nt, 8, 8, 8), dtype=np.complex128)

for conf in range(conf_count):
    confidx = nums[0] + nums[2] * conf

    gauge = io.readChromaQIOGauge(f'{folder}/{filename}/{filename}_{confidx}_coulomb_1e-9.scidac')
    print(f"gauge loaded... conf={confidx}")
    dslash.loadGauge(gauge)
    dslash0.loadGauge(gauge)
    fermion_delta_c1 = []
    fermion_delta_c2 = []
    fermion_delta_c3 = []
    dfermion_delta_c1 = []
    dfermion_delta_c2 = []
    dfermion_delta_c3 = []
    for delta in range(8):
        # inverse of D^+D
        propagator = invertStaggeredWithDelta(dslash, 0, delta)
        # apply dslash to each color
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
                    psi1 = fermion_delta_c[c1][bvd_idx]
                    psi2 = dfermion_delta_c[c1][bvd_idx]
                    phi1 = -psi2 + m * psi1
                    psi3 = fermion_delta_c[c1][B]
                    psi4 = dfermion_delta_c[c1][B]
                    phi2 = psi3 + m * psi4
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
    # Free propagators before next config
    np.save(f"{out_dir}/p_{conf}.npy", p_nt_a_b_d[conf, :, :, :, :])
    del fermion_delta_c, dfermion_delta_c
    print(f"  correlators conf={confidx} done")

# Project onto meson channels and save (same format as CorrelationP2P.py)
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
        np.save(f"{out_dir}/correlationp2p_{type}_{j}.npy", sum_res_lst)

correlation_func_combine = np.zeros(nt, dtype=np.complex128)
for channel in range(0, 20):
    correlation_func_combine = np.load(f"{out_dir}/correlationp2p_{channel}_0.npy")
    for i in range(1, len(all_signs_lst[channel])):
        correlation_func_combine = correlation_func_combine + np.load(f"{out_dir}/correlationp2p_{channel}_{i}.npy")
    np.save(f"{out_dir}/{basestr}_{channel}.npy", np.real(correlation_func_combine))

