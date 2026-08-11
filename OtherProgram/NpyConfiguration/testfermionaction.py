
"""
fermion action=22960.932550028654

"""
import numpy as np

white_noise = np.load("data/white_noise.npy")

print(np.shape(white_noise))
print(np.tensordot(np.conj(white_noise), white_noise, axes=5))