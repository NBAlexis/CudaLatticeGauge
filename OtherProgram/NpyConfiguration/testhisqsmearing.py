import numpy as np

from Converter.NpyConfigurationToCon import transformFromCon

epsilon = -0.1

clg_level1 = transformFromCon("data/cfg_100_clg_level1.con", 12, 12, 12, 12)
clg_level2 = transformFromCon("data/cfg_100_clg_level2.con", 12, 12, 12, 12)
clg_naik = transformFromCon("data/cfg_100_clg_naik.con", 12, 12, 12, 12)

pyquda_level1 = np.load("data/cfg_100_level1.npy")
pyquda_level2 = np.load("data/cfg_100_level2.npy")
pyquda_level2_with_eps = np.load("data/cfg_100_level2_eps.npy")
pyquda_naik_with_coefficient = np.load("data/cfg_100_naik.npy")
pyquda_naik_with_coefficient_and_eps = np.load("data/cfg_100_naik_eps.npy")

# npydata = transformFromCon("data/cfg_100_clg_naik.con", 12, 12, 12, 12)
# npydata = npydata * (-1/24)
# fatdata = np.load('data/cfg_100_naik.npy')

def staggeredPhase(x: int, y: int, z: int, t: int, l: int):
    site_list = np.array([x, y, z, t])
    return (-1)**(np.sum(site_list[:l]))

def staggeredPhaseMILC(x: int, y: int, z: int, t: int, l: int):
    site_list = np.array([t, x, y, z])
    l = (l + 1) % 4
    return (-1)**(np.sum(site_list[:l]))

def checkerror(d1, d2, withPhase = True):
    res = 0
    diff = 0
    shownexample = 10
    for x in range(12):
        for y in range(12):
            for z in range(12):
                for t in range(12):
                    for l in range(4):
                        m1 = d1[l, t, z, y, x]
                        m2 = d2[l, t, z, y, x] # * (staggeredPhaseMILC(x, y, z, t, l) if withPhase else 1)
                        m = m1 - m2
                        m = np.dot(m, m.transpose().conj())
                        deltam = m[0, 0] + m[1, 1] + m[2, 2]
                        if deltam > 1.0e-7:
                            diff = diff + 1
                            if shownexample > 0:
                                print(f"delta={deltam} ({x},{y},{z},{t})_{l},eta={staggeredPhaseMILC(x, y, z, t, l)}")
                                print(d1[l, t, z, y, x])
                                print(d2[l, t, z, y, x])
                                shownexample -= 1
                        res = res + m[0, 0] + m[1, 1] + m[2, 2]
    print(f"diff={diff}, ratio={diff/(12*12*12*12*4)}")
    return res

print(f"====== level1 ======== diff : {np.abs(checkerror(clg_level1, pyquda_level1))}")
print(f"====== level2 ======== diff : {np.abs(checkerror(clg_level2, pyquda_level2))}")
print(f"====== naik ======== diff : {np.abs(checkerror(clg_naik * (-1/24), pyquda_naik_with_coefficient))}")
print(f"====== level2 with eps ======== diff : {np.abs(checkerror(clg_level2 + clg_level1 * (epsilon/8), pyquda_level2_with_eps))}")
print(f"====== naik with eps ======== diff : {np.abs(checkerror(clg_naik * (-(1+epsilon)/24), pyquda_naik_with_coefficient_and_eps))}")
