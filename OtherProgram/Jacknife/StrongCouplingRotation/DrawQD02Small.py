import numpy as np

from Visualization import errorbar, LineStyle, MarkerStyle, LegendPosistion

data = np.loadtxt('Data/QD02Small.csv', delimiter=',')

print(np.shape(data))

errorbar([n * 0.0117 for n in range(0, 11)], [data[3,:]], [data[7,:]],
         xlabel="$a\\Omega$",
         ylabel="$P$",
         markers=[MarkerStyle.circle],
         linestyles=[LineStyle.dashed])

