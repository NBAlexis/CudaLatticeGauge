"""
for each sign_function, and delta

C(t) = sum_{x,y,A,B,c1,c2} s(A)s(B) S1_{x+A,y+B;c1,c2} S2*_{x+A+d,y+B+d;c1,c2}

so we pick up every same point from both propagators and do the sum

The quantity to be saved is p_{A,B,d}(t) = sum_{x,y}S1_{x+A,y+B;c1,c2} S2*_{x+A+d,y+B+d;c1,c2}

"""
import numpy as np

from MesonStructures import deltaidx_to_delta, delta_to_deltaidx

foldername = "G:\\mass\\c24p31\\"
m = 0.00944 * 2
nx = 24
nt = 48

for conf in range(51):
    conf_idx = 1000 + 20 * conf
    print("calculating conf = " + str(conf_idx))
    resarray = np.zeros((nt, 8, 8, 8), dtype=np.complex128)
    for A in range(8):
        for B in range(8):
            for delta in range(8):
                resarray[:, A, B, delta] = 0
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
                    psi1 = np.load(
                        foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{bvd_idx}_c{c1 + 1}_lexico.npy")
                    psi2 = np.load(
                        foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{bvd_idx}_dc{c1 + 1}_lexico.npy")
                    phi1 = -psi2 + m * psi1
                    psi3 = np.load(
                        foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{B}_c{c1 + 1}_lexico.npy")
                    psi4 = np.load(
                        foldername + f"propagator_s24t48_beta7.29_ml0.00944ms0.04721mc0.5555_{conf_idx}_delta{B}_dc{c1 + 1}_lexico.npy")
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
                    resarray[:, A, B, delta] = resarray[:, A, B, delta] + res
    # print(resarray)
    np.save(f"p_{conf_idx}.npy", resarray)
