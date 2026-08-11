import numpy as np

source_data = np.load("./data/Debug/source_delta3_t1_c1.npy")

all = 0
for x in range(24):
    for y in range(24):
        for z in range(24):
            for t in range(48):
                v = source_data[t, z, y, x, :]
                if np.sum(v * v) > 1.0e-6:
                    print(f"none-zero at xyzt = {x,y,z,t}: {v}")
                    all = all + 1
print(all)