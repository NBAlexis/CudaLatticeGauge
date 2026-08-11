import numpy as np

from Visualization import errorbar, LineStyle, MarkerStyle, LegendPosistion

data = np.loadtxt('Data/QD30.csv', delimiter=',')

nosite = 0
for i in range(14):
    for j in range(14):
        if (i - 7)** 2 + (j - 7)** 2 < 36:
            nosite += 1

print(nosite)

xilst = [0.25, 0.3, 0.35, 0.37, 0.4, 0.45, 0.5]
sep = len(xilst)

"""
errorbar([1/xilst[n] for n in range(sep)],
         [nosite*12*data[2*sep+0:2*sep+sep, 0], 14*14*12*data[7*sep+0:7*sep+sep, 0]],
         [nosite*12*data[3*sep+0:3*sep+sep, 0], 14*14*12*data[8*sep+0:8*sep+sep, 0]],
         xlabel="$1/\\xi$",
         ylabel="$\\chi _P$",
         linestyles=[LineStyle.dashdot, LineStyle.dashdot],
         markers=[MarkerStyle.circle, MarkerStyle.diamond],
         legends=["In", "All"])
errorbar([1/xilst[n] for n in range(sep)],
         [data[0+0:0+sep, 0], data[5*sep+0:5*sep+sep, 0]],
         [data[sep:sep+sep, 0], data[6*sep+0:6*sep+sep, 0]],
         xlabel="$1/\\xi$",
         ylabel="$P$",
         linestyles=[LineStyle.dashdot, LineStyle.dashdot],
         markers=[MarkerStyle.circle, MarkerStyle.diamond],
         legends=["In", "All"])
"""

# """
# errorbar([0.3 * n for n in range(11)], data[0+0:0+2,:], data[sep+0:sep+2,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.25$", "0.30"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.upperright)
# errorbar([0.3 * n for n in range(11)], data[5*sep+0:5*sep+2,:], data[6*sep+0:6*sep+2,:],
#          xlabel="$a\\Omega$",
#          ylabel="$P$",
#          legends=["$\\xi=0.25$", "0.30"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.upperright)
# errorbar([0.3 * n for n in range(11)], nosite*12*data[2*sep+0:2*sep+2,:], nosite*12*data[3*sep+0:3*sep+2,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.25$", "0.30"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# errorbar([0.3 * n for n in range(11)], 14*14*12*data[7*sep+0:7*sep+2,:], 14*14*12*data[8*sep+0:8*sep+2,:],
#          xlabel="$a\\Omega$",
#          ylabel="$\\chi_P$",
#          legends=["$\\xi=0.25$", "0.30"],
#          markers=[MarkerStyle.circle, MarkerStyle.diamond],
#          linestyles=[LineStyle.dashed, LineStyle.dashed],
#          legendpos=LegendPosistion.lowerright)
# """

"""
errorbar([0.3 * n for n in range(11)], data[0+2:0+7,:], data[sep+2:sep+7,:],
         xlabel="$a\\Omega$",
         ylabel="$P$",
         legends=["$\\xi=0.35$", "0.37", "0.40", "0.45", "0.50"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
         linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
         legendpos=LegendPosistion.lowerright)
errorbar([0.3 * n for n in range(11)], data[5*sep+2:5*sep+7,:], data[6*sep+2:6*sep+7,:],
         xlabel="$a\\Omega$",
         ylabel="$P$",
         legends=["$\\xi=0.35$", "0.37", "0.40", "0.45", "0.50"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
         linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
         legendpos=LegendPosistion.upperright)
errorbar([0.3 * n for n in range(11)], nosite*12*data[2*sep+2:2*sep+7,:], nosite*12*data[3*sep+2:3*sep+7,:],
         xlabel="$a\\Omega$",
         ylabel="$\\chi_P$",
         legends=["$\\xi=0.35$", "0.37", "0.40", "0.45", "0.50"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
         linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
         legendpos=LegendPosistion.lowerright)
errorbar([0.3 * n for n in range(11)], 14*14*12*data[7*sep+2:7*sep+7,:], 14*14*12*data[8*sep+2:8*sep+7,:],
         xlabel="$a\\Omega$",
         ylabel="$\\chi_P$",
         legends=["$\\xi=0.35$", "0.37", "0.40", "0.45", "0.50"],
         markers=[MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x],
         linestyles=[LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed, LineStyle.dashed],
         legendpos=LegendPosistion.lowerright)

"""

def PolyaFromTo(fromidx, toidx, showLgend):
    allmarkers = [MarkerStyle.circle, MarkerStyle.diamond, MarkerStyle.square, MarkerStyle.plus, MarkerStyle.x,
                  MarkerStyle.triangle_up, MarkerStyle.star, MarkerStyle.pentagon, MarkerStyle.triangle_down,
                  MarkerStyle.hexagon1, MarkerStyle.triangle_right, MarkerStyle.triangle_left, MarkerStyle.point, MarkerStyle.sp1,
                  MarkerStyle.sp2, MarkerStyle.sp3, MarkerStyle.sp4]
    markers = allmarkers[:(toidx - fromidx)]
    legends = [f"$\\xi={xilst[fromidx + n]}$" if 0 == n else f"{fromidx + n}" for n in range(toidx - fromidx)] if showLgend else None
    lines = [LineStyle.dashed for _ in range(toidx - fromidx)]
    errorbar([0.3 * n for n in range(11)], data[0*sep+fromidx:0*sep+toidx,:], data[1*sep+fromidx:1*sep+toidx,:],
             xlabel="$a\\Omega$",
             ylabel="$P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=LegendPosistion.upperright)
    errorbar([0.3 * n for n in range(11)], data[5*sep+fromidx:5*sep+toidx,:], data[6*sep+fromidx:6*sep+toidx,:],
             xlabel="$a\\Omega$",
             ylabel="$P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=LegendPosistion.upperright)
    errorbar([0.3 * n for n in range(11)], nosite*12*data[2*sep+fromidx:2*sep+toidx,:], nosite*12*data[3*sep+fromidx:3*sep+toidx,:],
             xlabel="$a\\Omega$",
             ylabel="$\\chi_P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=LegendPosistion.lowerright)
    errorbar([0.3 * n for n in range(11)], 14*14*12*data[7*sep+fromidx:7*sep+toidx,:], 14*14*12*data[8*sep+fromidx:8*sep+toidx,:],
             xlabel="$a\\Omega$",
             ylabel="$\\chi_P$",
             legends=legends,
             markers=markers,
             linestyles=lines,
             legendpos=LegendPosistion.lowerright)

PolyaFromTo(0, 7, False)