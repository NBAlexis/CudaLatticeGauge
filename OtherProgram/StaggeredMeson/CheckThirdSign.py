"""
for some channels, there are three different directions according to SO(3)
However, the third one has a different sign

But, why cannot I see it?
"""
import numpy as np
from matplotlib import pyplot as plt

pathstr = "./data/c32p31/correlationp2p"
checkchannels = [10, 11, 12, 13, 14, 15]

for i in range(len(checkchannels)):
    correlation_func1 = np.load(f"{pathstr}_{checkchannels[i]}_0.npy")
    correlation_func2 = np.load(f"{pathstr}_{checkchannels[i]}_1.npy")
    correlation_func3 = np.load(f"{pathstr}_{checkchannels[i]}_2.npy")
    correlation_func1ave = np.mean(np.real(correlation_func1), axis=0)
    correlation_func2ave = np.mean(np.real(correlation_func2), axis=0)
    correlation_func3ave = np.mean(np.real(correlation_func3), axis=0)
    plt.plot(correlation_func1ave)
    plt.plot(correlation_func2ave)
    plt.plot(correlation_func3ave)
    plt.show()
