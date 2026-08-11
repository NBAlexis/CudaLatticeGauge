import numpy as np

from Visualization import errorbar, LineStyle, MarkerStyle, LegendPosistion

data1 = np.loadtxt('Data/PolyaWD.csv', delimiter=',')
data2 = np.loadtxt('Data/ChiralWD.csv', delimiter=',')

print(np.shape(data1))
print(np.shape(data2))

# errorbar([n * 0.0117 for n in range(0, 11)],
#          [data1[0,:], data1[15,:]],
#          [data1[3,:], data1[18,:]],
#          xlabel="$a\\Omega$", ylabel="$P$",
#          legends=["In", "All"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed])
# errorbar([n * 0.0117 for n in range(0, 11)],
#          [data1[1,:], data1[16,:]],
#          [data1[4,:], data1[19,:]],
#          xlabel="$a\\Omega$", ylabel="$P$",
#          legends=["In", "All"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed])
# errorbar([n * 0.0117 for n in range(0, 11)],
#          [data1[2,:], data1[17,:]],
#          [data1[5,:], data1[20,:]],
#          xlabel="$a\\Omega$", ylabel="$P$",
#          legends=["In", "All"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed])


errorbar([n * 0.0117 for n in range(0, 11)],
         [data2[0,:], data2[2,:]],
         [data2[1,:], data2[3,:]],
         xlabel="$a\\Omega$", ylabel="$c$",
         legends=["In", "All"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond],
         linestyles=[LineStyle.dashed, LineStyle.dashed])

errorbar([n * 0.0117 for n in range(0, 11)],
         [data2[4,:], data2[6,:]],
         [data2[5,:], data2[7,:]],
         xlabel="$a\\Omega$", ylabel="$c$",
         legends=["In", "All"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond],
         linestyles=[LineStyle.dashed, LineStyle.dashed])
errorbar([n * 0.0117 for n in range(0, 11)],
         [data2[8,:], data2[10,:]],
         [data2[9,:], data2[11,:]],
         xlabel="$a\\Omega$", ylabel="$c$",
         legends=["In", "All"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond],
         linestyles=[LineStyle.dashed, LineStyle.dashed])