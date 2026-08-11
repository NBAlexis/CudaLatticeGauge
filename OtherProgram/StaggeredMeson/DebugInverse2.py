
"""

we have calculate:
phi1 = (D^+D)^-1 phi0
phi2 = D0 phi1
phi3 = D0 phi2

We expect:

phi2' = phi2 + m phi1 = (D^+)^-1 phi0
(-phi3 - m phi2 + m phi2 + m^2 phi1) = (-D0 (phi2 + m phi1) + m (phi2 + m phi1)) = phi0
"""
import numpy as np

phi1 = np.load("./data/Debug/inverse_1000_delta3_c1_lexico.npy")
# phi2 = np.load("./data/Debug/dinverse_1000_delta3_c1_lexico.npy")
phi3 = np.load("./data/Debug/ddinverse_1000_delta3_c1_lexico.npy")

m = 0.00944 * 2
phi0 = -phi3 + m * m * phi1

all = 0
for x in range(24):
    for y in range(24):
        for z in range(24):
            for t in range(48):
                v = phi0[t, z, y, x, :]
                if np.sum(v * v) > 1.0e-3:
                    print(f"none-zero at xyzt = {x,y,z,t}: {v}")
                    all = all + 1
print(all)