"""
to add more, see:
https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
https://matplotlib.org/stable/api/markers_api.html
"""
from enum import Enum

from matplotlib import pyplot as plt


class LineStyle(Enum):
    none = 0
    solid = 1
    dotted = 2
    dashed = 3
    dashdot = 4
    dashdotdot = 5
    loosely_dotted = 6
    densely_dotted = 7
    loosely_dashed = 8
    densely_dashed = 9
    loosely_dashdotted = 10
    densely_dashdotted = 11
    loosely_dashdotdotted = 12
    densely_dashdotdotted = 13


class MarkerStyle(Enum):
    none = 0
    point = 1
    circle = 2
    square = 3
    plus = 4
    x = 5
    diamond = 6
    triangle_down = 7
    triangle_up = 8
    star = 9
    hexagon1 = 10
    triangle_left = 11
    triangle_right = 12
    pentagon = 13
    sp1 = 14
    sp2 = 15
    sp3 = 16
    sp4 = 17

class LegendPosistion(Enum):
    best = 0
    upperright = 1
    upperleft = 2
    lowerleft = 3
    lowerright = 4
    right = 5
    centerleft = 6
    centerright = 7
    lowercenter = 8
    uppercenter = 9
    center = 10


linestyle_dic = {
    LineStyle.none: 'none',
    LineStyle.solid: 'solid',
    LineStyle.dotted: 'dotted',
    LineStyle.dashed: 'dashed',
    LineStyle.dashdot: 'dashdot',
    LineStyle.dashdotdot: (0, (3, 5, 1, 5, 1, 5)),
    LineStyle.loosely_dotted: (0, (1, 10)),
    LineStyle.densely_dotted: (0, (1, 1)),
    LineStyle.loosely_dashed: (0, (5, 10)),
    LineStyle.densely_dashed: (0, (5, 1)),
    LineStyle.loosely_dashdotted: (0, (3, 10, 1, 10)),
    LineStyle.densely_dashdotted: (0, (3, 1, 1, 1)),
    LineStyle.loosely_dashdotdotted: (0, (3, 10, 1, 10, 1, 10)),
    LineStyle.densely_dashdotdotted: (0, (3, 1, 1, 1, 1, 1)),
}

marker_dic = {
    MarkerStyle.none: ' ',
    MarkerStyle.point: '.',
    MarkerStyle.circle: 'o',
    MarkerStyle.square: 's',
    MarkerStyle.plus: '+',
    MarkerStyle.x: 'x',
    MarkerStyle.diamond: 'd',
    MarkerStyle.triangle_down: 'v',
    MarkerStyle.triangle_up: '^',
    MarkerStyle.star: '*',
    MarkerStyle.hexagon1: 'h',
    MarkerStyle.triangle_left: '<',
    MarkerStyle.triangle_right: '>',
    MarkerStyle.pentagon: 'p',
    MarkerStyle.sp1: '1',
    MarkerStyle.sp2: '2',
    MarkerStyle.sp3: '3',
    MarkerStyle.sp4: '4',
}

legendpos_dic = {
    LegendPosistion.best: 'best',
    LegendPosistion.upperright: 'upper right',
    LegendPosistion.upperleft: 'upper left',
    LegendPosistion.lowerleft: 'lower left',
    LegendPosistion.lowerright: 'lower right',
    LegendPosistion.right: 'right',
    LegendPosistion.centerleft: 'center left',
    LegendPosistion.centerright: 'center right',
    LegendPosistion.lowercenter: 'lower center',
    LegendPosistion.uppercenter: 'upper center',
    LegendPosistion.center: 'center',
}


def errorbar(x, ylst, errlst, **kwargs):
    xlabel = kwargs.pop('xlabel', None)
    ylabel = kwargs.pop('ylabel', None)
    ylog = kwargs.pop('ylog', False)
    markers = kwargs.pop('markers', None)
    linestyles = kwargs.pop('linestyles', None)
    markerSize = kwargs.pop('markerSize', 5.0)
    lineWidth = kwargs.pop('lineWidth', 1.0)
    legends = kwargs.pop('legends', None)
    legendpos = kwargs.pop('legendpos', LegendPosistion.best)
    savefile = kwargs.pop('savefile', None)
    fig = plt.figure(figsize=(5, 4), layout="constrained")
    axis = fig.add_subplot(111)
    for i in range(len(ylst)):
        marker = None
        linestyle = None
        showlegends = legends[i] if legends is not None and legends[i] is not None else None
        if markers is not None:
            marker = marker_dic[markers[i]]
        if linestyles is not None:
            linestyle = linestyle_dic[linestyles[i]]
        if marker is None:
            if linestyle is None:
                axis.errorbar(x, ylst[i], yerr=errlst[i], label=showlegends, linewidth=lineWidth)
            else:
                axis.errorbar(x, ylst[i], yerr=errlst[i], linestyle=linestyle, label=showlegends, linewidth=lineWidth)
        else:
            if linestyle is None:
                axis.errorbar(x, ylst[i], yerr=errlst[i], marker=marker, markersize=markerSize, label=showlegends, linewidth=lineWidth)
            else:
                axis.errorbar(x, ylst[i], yerr=errlst[i], marker=marker, markersize=markerSize, linestyle=linestyle, label=showlegends, linewidth=lineWidth)
    if xlabel is not None:
        axis.set_xlabel(xlabel)
    if ylabel is not None:
        axis.set_ylabel(ylabel)
    if ylog:
        axis.set_yscale('log')
    if legends is not None:
        axis.legend(loc=legendpos_dic[legendpos])
    if savefile is None:
        plt.show()
    else:
        plt.savefig(savefile)

def history(yarrays, **kwargs):
    showrange = kwargs.pop('range', None)
    legends = kwargs.pop('legends', None)
    xstart = kwargs.pop('xstart', 0)
    xlabel = kwargs.pop('xlabel', None)
    ylabel = kwargs.pop('ylabel', None)
    fig = plt.figure(figsize=(10, 6))
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(1, 2, width_ratios=(4, 1),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    # Create the Axes.
    ax = fig.add_subplot(gs[0, 0])
    ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
    xpoints = [i + xstart + 1 for i in range(len(yarrays[0]))]
    for j in range(len(yarrays)):
        showlegends = legends[j] if legends is not None and legends[j] is not None else None
        ax.scatter(xpoints, yarrays[j], s=1, label=showlegends)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if legends is not None:
        ax.legend()

    miny = min([min(y) for y in yarrays])
    maxy = max([max(y) for y in yarrays])
    if showrange is not None:
        miny = showrange[0]
        maxy = showrange[1]
    for j in range(len(yarrays)):
        ax_histy.hist(yarrays[j], bins=50, range=[miny, maxy], histtype='step', orientation='horizontal')
    plt.show()
