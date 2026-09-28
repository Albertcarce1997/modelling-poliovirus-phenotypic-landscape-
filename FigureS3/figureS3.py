"""
Figure S3 - how far the haplotypes lie from the most abundant one.

Each nucleotide haplotype is placed in a bin by the number of nucleotide
differences separating it from HT-1, the most abundant haplotype.  Bars give
the number of haplotypes in each bin and carry that count above them; the line
gives the share of all reads those haplotypes represent, on the right axis.
HT-1 itself is not drawn but stays in the read denominator.

Run:      python figureS3.py
Data:     haplotype_mutation_counts.csv
Output:   figureS3.pdf, .svg and .png, and figureS3_values.csv
Requires: python >= 3.9, matplotlib, pandas
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))

DATA = 'haplotype_mutation_counts.csv'
STEM = 'figureS3'
TABLE = 'figureS3_values.csv'
FORMATS = ('pdf', 'svg', 'png')
PNG_DPI = 300
FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']

DROP_ZERO = True                   # HT-1 itself, at distance 0
C_BAR = '#2f6f73'                  # haplotype count
C_LINE = '#7a4f9a'                 # share of reads
C_DOT = '#111111'
C_GRID = '#e6e6e6'
C_AXIS = '#222222'
C_TEXT = '#333333'

BAR_WIDTH = 0.55
LINE_WIDTH, DOT_SIZE = 1.6, 22.0
YLIM_LEFT, YSTEP_LEFT = (0, 800), 200
YLIM_RIGHT, YSTEP_RIGHT = (0, 30), 10
FIG_W, FIG_H = 8.3, 4.2
FS_TICK, FS_LAB, FS_ANN = 9.0, 10.5, 8.5

XLABEL = 'Number of mutations'
YLABEL_LEFT = 'Number of haplotypes'
YLABEL_RIGHT = 'Percent of total reads'
LABEL_BAR = 'Haplotype count'
LABEL_LINE = 'Percent of total reads'

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': FONT,
                     'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 1.0})


def main():
    path = os.path.join(HERE, DATA)
    if not os.path.exists(path):
        sys.exit('%s not found next to the script' % DATA)
    d = pd.read_csv(path)
    total_reads = d['Reads'].sum()
    g = (d.groupby('Mutations_vs_HT1')
           .agg(Haplotypes=('Haplotype', 'size'), Reads=('Reads', 'sum'))
           .reset_index())
    g['Percent_of_reads'] = 100.0 * g['Reads'] / total_reads
    if DROP_ZERO:
        g = g[g['Mutations_vs_HT1'] > 0]
    g = g.sort_values('Mutations_vs_HT1').reset_index(drop=True)
    g.round(4).to_csv(os.path.join(HERE, TABLE), index=False)

    x = np.arange(len(g))
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for t in range(YLIM_LEFT[0], YLIM_LEFT[1] + 1, YSTEP_LEFT):
        ax.axhline(t, color=C_GRID, lw=1.0, zorder=0)

    ax.bar(x, g['Haplotypes'], color=C_BAR, width=BAR_WIDTH, zorder=2)
    for xi, n in zip(x, g['Haplotypes']):
        ax.annotate('%d' % n, (xi, n), textcoords='offset points', xytext=(0, 4),
                    ha='center', fontsize=FS_ANN, color=C_AXIS)

    ax2 = ax.twinx()
    ax2.plot(x, g['Percent_of_reads'], color=C_LINE, lw=LINE_WIDTH, zorder=3)
    ax2.scatter(x, g['Percent_of_reads'], s=DOT_SIZE, c=C_DOT, zorder=4)

    ax.set_xlim(-0.6, len(g) - 0.4)
    ax.set_ylim(*YLIM_LEFT)
    ax2.set_ylim(*YLIM_RIGHT)
    ax.set_xticks(x)
    ax.set_xticklabels(g['Mutations_vs_HT1'])
    ax.set_yticks(range(YLIM_LEFT[0], YLIM_LEFT[1] + 1, YSTEP_LEFT))
    ax2.set_yticks(range(YLIM_RIGHT[0], YLIM_RIGHT[1] + 1, YSTEP_RIGHT))
    ax2.set_yticklabels(['%d%%' % t for t in
                         range(YLIM_RIGHT[0], YLIM_RIGHT[1] + 1, YSTEP_RIGHT)])
    ax.set_xlabel(XLABEL, fontsize=FS_LAB, color=C_AXIS)
    ax.set_ylabel(YLABEL_LEFT, fontsize=FS_LAB, color=C_AXIS)
    ax2.set_ylabel(YLABEL_RIGHT, fontsize=FS_LAB, color=C_AXIS)
    for a in (ax, ax2):
        a.tick_params(labelsize=FS_TICK, width=1.0, length=0, colors=C_TEXT)
    for sp in ('top',):
        ax.spines[sp].set_visible(False)
        ax2.spines[sp].set_visible(False)
    for sp in ('left', 'bottom', 'right'):
        ax.spines[sp].set_color(C_AXIS)
        ax2.spines[sp].set_color(C_AXIS)

    ax.legend(handles=[Patch(facecolor=C_BAR, edgecolor='none', label=LABEL_BAR),
                       Line2D([], [], color=C_LINE, lw=LINE_WIDTH, marker='o',
                              markerfacecolor=C_DOT, markeredgecolor=C_DOT,
                              markersize=4, label=LABEL_LINE)],
              loc='upper right', frameon=False, fontsize=FS_ANN, ncol=2)

    fig.tight_layout()
    for f in FORMATS:
        fig.savefig(os.path.join(HERE, '%s.%s' % (STEM, f)), format=f,
                    dpi=PNG_DPI if f == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    print(g.to_string(index=False))
    print('\nwritten to %s.{%s} and %s' % (STEM, ','.join(FORMATS), TABLE))


if __name__ == '__main__':
    main()
