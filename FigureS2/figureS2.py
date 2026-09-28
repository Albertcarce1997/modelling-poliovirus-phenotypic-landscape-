"""
Figure S2 - mean quality score along the VP1 gene.

One point per VP1 nucleotide position, giving the mean Phred quality score of
the base calls covering it across every retained amino acid haplotype.

Run:      python figureS2.py
Data:     position_qscores.csv
Output:   figureS2.pdf, .svg and .png
Requires: python >= 3.9, matplotlib, pandas
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

DATA = 'position_qscores.csv'
STEM = 'figureS2'
FORMATS = ('pdf', 'svg', 'png')
PNG_DPI = 300
FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']

C_LINE = '#2f6f73'                 # the mean trace
C_GRID = '#e6e6e6'
C_AXIS = '#222222'
C_TEXT = '#333333'
LINE_WIDTH = 1.4

YLIM, YSTEP = (0, 50), 10
N_XTICKS = 7                       # evenly spaced, first and last included
FIG_W, FIG_H = 8.3, 3.4
FS_TICK, FS_LAB = 9.0, 10.5

XLABEL = 'VP1 nucleotide position in read'
YLABEL = 'Qscore'

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': FONT,
                     'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 1.0})


def main():
    path = os.path.join(HERE, DATA)
    if not os.path.exists(path):
        sys.exit('%s not found next to the script' % DATA)
    d = pd.read_csv(path).sort_values('vp1_nt_position')
    x = d['vp1_nt_position'].to_numpy()
    y = d['mean_qscore'].to_numpy()

    ticks = np.unique(np.rint(np.linspace(x.min(), x.max(), N_XTICKS)).astype(int))
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for t in ticks:
        ax.axvline(t, color=C_GRID, lw=1.0, zorder=0)
    for t in range(YLIM[0], YLIM[1] + 1, YSTEP):
        ax.axhline(t, color=C_GRID, lw=1.0, zorder=0)

    ax.plot(x, y, color=C_LINE, lw=LINE_WIDTH, solid_joinstyle='round',
            solid_capstyle='round', zorder=3)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(*YLIM)
    ax.set_xticks(ticks)
    ax.set_yticks(range(YLIM[0], YLIM[1] + 1, YSTEP))
    ax.set_xlabel(XLABEL, fontsize=FS_LAB, color=C_AXIS)
    ax.set_ylabel(YLABEL, fontsize=FS_LAB, color=C_AXIS)
    ax.tick_params(labelsize=FS_TICK, width=1.0, length=0, colors=C_TEXT)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(C_AXIS)

    fig.tight_layout()
    for f in FORMATS:
        fig.savefig(os.path.join(HERE, '%s.%s' % (STEM, f)), format=f,
                    dpi=PNG_DPI if f == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print('positions %d   overall mean Q %.2f   lowest %.2f at position %d'
          % (len(d), y.mean(), y.min(), int(x[y.argmin()])))
    print('written to %s.{%s}' % (STEM, ','.join(FORMATS)))


if __name__ == '__main__':
    main()
