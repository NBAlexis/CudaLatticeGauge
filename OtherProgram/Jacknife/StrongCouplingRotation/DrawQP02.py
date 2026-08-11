import numpy as np
from matplotlib import pyplot as plt

from Visualization import errorbar, LineStyle, MarkerStyle, LegendPosistion

data = np.loadtxt('Data/QP02.csv', delimiter=',')

nosite = 0
for i in range(12):
    for j in range(12):
        if (i - 6)** 2 + (j - 6)** 2 < 36:
            nosite += 1

print(nosite)
xilst = [0.0227, 0.02275, 0.0228, 0.02285, 0.0229, 0.02295, 0.0230, 0.02305, 0.0231, 0.02315, 0.0232, 0.02325, 0.0233, 0.02335, 0.0234]
sep = len(xilst)
print(np.shape(data))

# errorbar([1/xilst[n] for n in range(sep)],
#          [12*12*12*data[7*sep+0:7*sep+sep, 0]],
#          [12*12*12*data[8*sep+0:8*sep+sep, 0]],
#          xlabel="$1/\\xi$",
#          ylabel="$\\chi _P$",
#          linestyles=[LineStyle.dashed],
#          markers=[MarkerStyle.circle],
#          savefile="Data/Fig/QP02_PolyakovChi.pdf")
# errorbar([1/xilst[n] for n in range(sep)],
#          [data[5*sep+0:5*sep+sep, 0]],
#          [data[6*sep+0:6*sep+sep, 0]],
#          xlabel="$1/\\xi$",
#          ylabel="$P$",
#          linestyles=[LineStyle.dashed],
#          markers=[MarkerStyle.circle],
#          savefile="Data/Fig/QP02_Polyakov.pdf")


def ExportPolyaAndChi(idx, left, right, filename):
    errorbar([left + n for n in range(right - left)], [data[5*sep+idx,left:right]], [data[6*sep+idx,left:right]],
             xlabel="$a\\Omega$",
             ylabel="$P$",
             linestyles=[LineStyle.dashed],
             markers=[MarkerStyle.circle],
             savefile=f"Data/Fig/QP02_{filename}_Polyakov.pdf")
    errorbar([left + n for n in range(right - left)], [12*12*12*data[7*sep+idx,left:right]], [12*12*12*data[8*sep+idx,left:right]],
             xlabel="$a\\Omega$",
             ylabel="$\\chi_P$",
             linestyles=[LineStyle.dashed],
             markers=[MarkerStyle.circle],
             savefile=f"Data/Fig/QP02_{filename}_PolyakovChi.pdf")

def PolyaFromTo(fromidx, toidx, left, right, showLgend, polyapos, susppos):
    allmarkers = [MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x,
                  MarkerStyle.triangle_up, MarkerStyle.star, MarkerStyle.pentagon, MarkerStyle.triangle_down,
                  MarkerStyle.hexagon1, MarkerStyle.triangle_right, MarkerStyle.triangle_left, MarkerStyle.point, MarkerStyle.sp1,
                  MarkerStyle.sp2, MarkerStyle.sp3, MarkerStyle.sp4]
    markers = allmarkers[:(toidx - fromidx)]
    legends = [f"$\\xi={xilst[fromidx + n]}$" if 0 == n else f"{xilst[fromidx + n]}" for n in range(toidx - fromidx)] if showLgend else None
    lines = [LineStyle.dashed for _ in range(toidx - fromidx)]
    errorbar([left + n for n in range(right - left)], data[0*sep+fromidx:0*sep+toidx,left:right], data[1*sep+fromidx:1*sep+toidx,left:right],
             xlabel="$a\\Omega$",
             ylabel="$P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=polyapos)
    errorbar([left + n for n in range(right - left)], data[5*sep+fromidx:5*sep+toidx,left:right], data[6*sep+fromidx:6*sep+toidx,left:right],
             xlabel="$a\\Omega$",
             ylabel="$P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=polyapos)
    errorbar([left + n for n in range(right - left)], nosite*12*data[2*sep+fromidx:2*sep+toidx,left:right], nosite*12*data[3*sep+fromidx:3*sep+toidx,left:right],
             xlabel="$a\\Omega$",
             ylabel="$\\chi_P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=susppos)
    errorbar([left + n for n in range(right - left)], 14*14*12*data[7*sep+fromidx:7*sep+toidx,left:right], 14*14*12*data[8*sep+fromidx:8*sep+toidx,left:right],
             xlabel="$a\\Omega$",
             ylabel="$\\chi_P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=susppos)

xilst = [2270, 2275, 2280, 2285, 2290, 2295, 2300, 2305, 2310, 2315, 2320, 2325, 2330, 2335, 2340]
for i in range(len(xilst)):
    ExportPolyaAndChi(i, 0, 31, xilst[i])

# PolyaFromTo(13, 16, 0, 5, False, LegendPosistion.upperright, LegendPosistion.upperright)

# errorbar(range(11), data[0+4:0+6,:], data[11+4:11+6,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.022$", "0.024"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), data[55+4:55+6,:], data[66+4:66+6,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.022$", "0.024"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), data[0+6:0+11,:], data[11+6:11+11,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.026", "0.028", "0.030", "0.032", "0.034"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), data[55+6:55+11,:], data[66+6:66+11,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.026", "0.028", "0.030", "0.032", "0.034"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)

# errorbar(range(11), nosite*12*data[22+0:22+4,:], nosite*12*data[33+0:33+4,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.014$", "0.016", "0.018", "0.020"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), 14*14*12*data[77+0:77+4,:], 14*14*12*data[88+0:88+4,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.014$", "0.016", "0.018", "0.020"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), nosite*12*data[22+4:22+6,:], nosite*12*data[33+4:33+6,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.022$", "0.024"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.upperright)
# errorbar(range(11), 14*14*12*data[77+4:77+6,:], 14*14*12*data[88+4:88+6,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.022$", "0.024"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.upperright)
# errorbar(range(11), nosite*12*data[22+6:22+11,:], nosite*12*data[33+6:33+11,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.026", "0.028", "0.030", "0.032", "0.034"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar(range(11), 14*14*12*data[77+6:77+11,:], 14*14*12*data[88+6:88+11,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.026", "0.028", "0.030", "0.032", "0.034"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
#          linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)

# errorbar(range(11), [nosite*12*data[22+4,:], 14*14*12*data[77+4,:]],
#                    [nosite*12*data[33+4,:], 14*14*12*data[88+4,:]],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["In", "All"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.upperright)
