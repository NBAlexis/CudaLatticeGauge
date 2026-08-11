import numpy as np

from Converter.NpyConfigurationToCon import transformStaggeredFermionFromCon

prop = np.load("data/pointsource_propagator_cfg_100.npy")
clg_prop_r = transformStaggeredFermionFromCon("data/clg_pointsource_cfg_100_r_double.con", 12, 12, 12, 12, True)
clg_prop_g = transformStaggeredFermionFromCon("data/clg_pointsource_cfg_100_g_double.con", 12, 12, 12, 12, True)
clg_prop_b = transformStaggeredFermionFromCon("data/clg_pointsource_cfg_100_b_double.con", 12, 12, 12, 12, True)

def checkerror(d1, d2_r, d2_g, d2_b):
    res = 0
    diff = 0
    shownexample = 10
    for x in range(12):
        for y in range(12):
            for z in range(12):
                for t in range(12):
                    m1 = d1[t, z, y, x]
                    m2 = np.array([d2_r[t, z, y, x], d2_g[t, z, y, x], d2_b[t, z, y, x]])
                    m2 = np.transpose(m2)
                    m = m1 - m2
                    m = np.dot(m, m.transpose().conj())
                    deltam = m[0, 0] + m[1, 1] + m[2, 2]
                    if deltam > 1.0e-7:
                        diff = diff + 1
                        if shownexample > 0:
                            print(f"delta={deltam} ({x},{y},{z},{t})")
                            print(m1)
                            print(m2)
                            shownexample -= 1
                    res = res + m[0, 0] + m[1, 1] + m[2, 2]
    print(f"diff={diff}, ratio={diff/(12*12*12*12)}")
    return res

print(checkerror(prop, clg_prop_r, clg_prop_g, clg_prop_b))
