import numpy as np

dslashres = np.load("data/dslash_identity_point.npy")
lx, ly, lz, lt = 12, 12, 12, 12

nonezero = 0
for x in range(lx):
    for y in range(ly):
        for z in range(lz):
            for t in range(lt):
                if (np.abs(dslashres[t,z,y,x,0])**2) + (np.abs(dslashres[t,z,y,x,1])**2) + (np.abs(dslashres[t,z,y,x,2])**2) > 1.0e-15:
                    print(f"non-zero at {x},{y},{z},{t} = {dslashres[t,z,y,x]}")
                    nonezero = nonezero + 1

print(nonezero)
