"""
we have
psi1 = (D^+D)^-1 source
and
psi2 = D0 psi1
we want:

phi1 = D^+ psi1 = (-D0 + m) psi1 = -psi2 + m psi1
phi2 = D   psi1 = ( D0 + m) psi1 =  psi2 + m psi1

and then we need sink, which is equivalent to do sum directly

u,v [t, sink_v, source_v, sink_c, source_c]
"""
import numpy as np

foldername = "G:\\mass\\c24p31\\"
m = 0.00944 * 2
nx = 24
nt = 48

resarray1 = np.zeros((nt, 8, 8, 3, 3), dtype=np.complex128)
resarray2 = np.zeros((nt, 8, 8, 3, 3), dtype=np.complex128)

for conf in range(51):
    conf_idx = 1000 + 20 * conf
    print("calculating conf = " + str(conf_idx))
    for delta1 in range(8):
        for c1 in range(3):
            psi1 = np.load(foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{delta1}_c{c1 + 1}_lexico.npy")
            psi2 = np.load(foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{delta1}_dc{c1 + 1}_lexico.npy")
            phi1 = -psi2 + m * psi1
            phi2 =  psi2 + m * psi1
            for delta2 in range(8):
                for c2 in range(3):
                    """
                    B is delta1
                    A is delta2
                    c2 is c1
                    c1 is c2
                    """
                    for t in range(nt):
                        # t,z,y,x
                        # take A,c1 which is delta2,c2
                        sum1array = phi1[t, ((delta2>>2)&1):nx:2, ((delta2>>1)&1):nx:2, (delta2&1):nx:2, c2]
                        sum2array = phi2[t, ((delta2>>2)&1):nx:2, ((delta2>>1)&1):nx:2, (delta2&1):nx:2, c2]
                        resarray1[t, delta2, delta1, c2, c1] = np.sum(sum1array)
                        resarray2[t, delta2, delta1, c2, c1] = np.sum(sum2array)
    np.save(f"data/u{conf_idx}.npy", resarray1)
    np.save(f"data/v{conf_idx}.npy", resarray2)




