import numpy as np

from Converter.NpyConfigurationToCon import transformStaggeredFermionFromCon

dslashres = np.load("data/dslash_cfg_100_point.npy")
npres = transformStaggeredFermionFromCon("data/clg_d_even_on_point_cfg_100_r_double.con", 12, 12, 12, 12, True)
print(np.shape(npres))
lx, ly, lz, lt = 12, 12, 12, 12

def checkerrorvector(d1, d2):
    res = 0
    diff = 0
    shownexample = 10
    for x in range(12):
        for y in range(12):
            for z in range(12):
                for t in range(12):
                    for c in range(3):
                        m1 = d1[t, z, y, x, c]
                        m2 = d2[t, z, y, x, c]
                        deltam = np.abs(m1 - m2)
                        if deltam > 1.0e-7:
                            diff = diff + 1
                            if shownexample > 0:
                                print(f"delta={deltam} ({x},{y},{z},{t})_{c}")
                                print(d1[t, z, y, x, c])
                                print(d2[t, z, y, x, c])
                                shownexample -= 1
                        res = res + deltam
    print(f"diff={diff}, ratio={diff/(12*12*12*12*3)}")
    return res

print(checkerrorvector(dslashres, npres))

