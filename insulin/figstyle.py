"""Shared figure style: no in-figure titles (the caption carries the bold title),
lowercase a)/b)/c) panel labels, and consistent fonts/spines."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def setup():
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.dpi': 140,
                         'axes.titlesize': 10, 'legend.frameon': False})


def panel(ax, letter, subtitle=''):
    """Panel label 'a) short subtitle' at the top-left of the axes, in the
    (otherwise unused) title slot so it never collides with data. The subtitle
    is a brief complement to the caption, not a repeat of it."""
    txt = ('%s) %s' % (letter, subtitle)).rstrip()
    ax.set_title(txt, loc='left', fontweight='bold', fontsize=9, pad=4)
