import numpy as np

from MesonStructures import all_signs
all_signs_lst = all_signs()

pathstrlst = ["c24p22", "c24p31", "c24p31s", "c32p22", "c32p31", "e32p31", "g32p32"]

for pathstr in pathstrlst:
    for channel in range(0, 20):
        correlation_func = np.load(f"data/{pathstr}/correlationp2p_{channel}_0.npy")
        for i in range(1, len(all_signs_lst[channel])):
            correlation_func = correlation_func + np.load(f"data/{pathstr}/correlationp2p_{channel}_{i}.npy")
        np.save(f"data/correlators/{pathstr}/{pathstr}_{channel}.npy", np.real(correlation_func))
    print(f"{pathstr} done")

