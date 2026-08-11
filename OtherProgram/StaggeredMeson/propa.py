import numpy as np
import sys

nx = 32
nt = 64
mass = 0.00743
# mass = 0.00579

args = sys.argv[1:]
nums = [int(arg) for arg in args]

print(f"start_idx = {nums[0]} count = {nums[1]} stride = {nums[2]}")

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

init(
    grid_size=None,
    latt_size=[nx, nx, nx, nt],
    resource_path=".cache"
)
latt_info = core.LatticeInfo([nx, nx, nx, nt], -1, 1.0)

dslash = core.getHISQ(latt_info, mass, 1e-12, 5000, 0.0)
dslash.invert_param.solution_type = QudaSolutionType.QUDA_MATDAG_MAT_SOLUTION
dslash.invert_param.solve_type = QudaSolveType.QUDA_NORMOP_SOLVE
# dslash.invert_param.inv_type = QudaInverterType.QUDA_BICGSTAB_INVERTER

#zero mass dslash to do D0 on (D^+D)^-1 phi
dslash0 = core.getHISQ(latt_info, 0.0, 1e-12, 5000, 0.0)

for conf in range(nums[1]):
    confidx = nums[0] + nums[2] * conf
    
    gauge = io.readChromaQIOGauge(f'/public/share/ybyang/caimc/CloveronHISQ/s32t64_beta7.54_ml0.00743ms0.03715mc0.4371/s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_coulomb_1e-9.scidac')
    # gauge = io.readChromaQIOGauge(f'/public/share/ybyang/caimc/CloveronHISQ/s32t64_beta7.75_ml0.00579ms0.02895mc0.34/s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_coulomb_1e-9.scidac')
    print(f"gauge loaded...conf={confidx}")
    dslash.loadGauge(gauge)
    dslash0.loadGauge(gauge)
    for delta in range(8):
        # inverse of D^+D
        propagator = invertStaggeredWithDelta(dslash, 0, delta)
        # apply dslash to each color
        fermion_0 = propagator.getFermion(0)
        fermion_1 = propagator.getFermion(1)
        fermion_2 = propagator.getFermion(2)
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_c1_lexico.npy", fermion_0.lexico())
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_c2_lexico.npy", fermion_1.lexico())
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_c3_lexico.npy", fermion_2.lexico())        
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_c1_lexico.npy", fermion_0.lexico())
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_c2_lexico.npy", fermion_1.lexico())
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_c3_lexico.npy", fermion_2.lexico())        
        dfermion_0 = dslash0.dslash(fermion_0)
        dfermion_1 = dslash0.dslash(fermion_1)
        dfermion_2 = dslash0.dslash(fermion_2)
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_dc1_lexico.npy", dfermion_0.lexico())
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_dc2_lexico.npy", dfermion_1.lexico())
        np.save(f"./data/propagator_s32t64_beta7.54_ml0.00743ms0.03715mc0.4371_{confidx}_delta{delta}_dc3_lexico.npy", dfermion_2.lexico())
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_dc1_lexico.npy", dfermion_0.lexico())
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_dc2_lexico.npy", dfermion_1.lexico())
        # np.save(f"./data/propagator_s32t64_beta7.75_ml0.00579ms0.02895mc0.34_{confidx}_delta{delta}_dc3_lexico.npy", dfermion_2.lexico())        
        print(f"result conf={confidx}-delta={delta} saved")
