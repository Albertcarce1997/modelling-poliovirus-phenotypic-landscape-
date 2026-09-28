"""
Figure 2, panel E - receptor against pocket-factor binding change.

Every Barcelona structural haplotype is one dot: its ddG against the CD155 receptor on
the horizontal axis and against the pocket factor on the vertical axis.  Dashed
lines mark Sabin 2 at zero, dotted lines the +/- 1 REU criterion.  Sets that
cross the criterion for both interactors are drawn on top of the cloud, coloured
where the pocket-factor change also exceeds 2 REU and black otherwise, so that
they arrive as their own groups in the vector output and can be labelled by hand.

Run:      python figure2_panelE.py
Data:     VDPV2_VP1_binding_energy_raw_data.xlsx  (or $VDPV2_DATA)
Output:   figures/, as PDF, SVG and PNG with text kept as text
Requires: python >= 3.9, matplotlib >= 3.6, numpy, pandas, openpyxl
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ###########################################################################
# CONFIG
# ###########################################################################

OUTDIR  = 'figures'
FORMATS = ('pdf', 'svg', 'png')     # png is a preview only
PNG_DPI = 200

FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']

DATASET = 'Barcelona sewage (haplotype-resolved)'
X_INTERACTOR = 'CD155 receptor'
Y_INTERACTOR = 'Pocket factor (palmitate)'

THRESHOLD = 1.0                     # REU; dotted lines on both axes
STRONG    = 2.0                     # pocket-factor ddG above which a crossing
                                    # set is coloured rather than black

C_CLOUD = '#C8C8C8'                 # every structural haplotype
C_BLACK = '#1A1A1A'                 # crosses the criterion for both interactors
C_AXIS  = '#333333'
C_GRID  = '#D9D9D9'
C_ZERO  = '#333333'

# one colour per highlighted set, keyed by the haplotype it is named after
COLOURS = {'HT-837': '#F4A0A8', 'HT-485': '#8FD9B6', 'HT-877': '#A8CC4A',
           'HT-442': '#E87DC4', 'HT-294': '#6FB6E8', 'HT-886': '#9B6FC4',
           'HT-385': '#3B7FD4', 'HT-399': '#C8C800'}

DOT_CLOUD, DOT_MARK = 26.0, 34.0    # marker area, pt^2
ALPHA_CLOUD = 0.75
EDGE_MARK = 0.0                     # outline width on the marked dots

XLIM, YLIM = (-3.8, 3.8), (-1.8, 11.2)
XTICKS, YTICKS = (-2, 0, 2), (0, 4, 8)

FIG_W, FIG_H = 5.2, 3.9             # inches
FS_TICK, FS_LAB = 8.0, 9.0

XLABEL = 'Receptor ΔΔG (REU)'
YLABEL = 'Pocket factor ΔΔG (REU)'

WRITE_POSITIONS = True              # csv of every crossing set, for labelling

# ###########################################################################
# end of CONFIG
# ###########################################################################

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': FONT,
                     'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 0.6})


HERE = os.path.dirname(os.path.abspath(__file__))


def _find_data():
    return os.environ.get('VDPV2_DATA') or os.path.join(HERE, 'VDPV2_VP1_binding_energy_raw_data.xlsx')


def load():
    d = pd.read_excel(_find_data(), sheet_name='Data')
    d = d[(d['Dataset'] == DATASET) & (d['Structural_haplotype'] != 'Sabin 2 (reference)')]
    x = d[d['Interactor'] == X_INTERACTOR][['Structural_haplotype', 'sHT', 'Haplotypes', 'ddG_REU']]
    y = d[d['Interactor'] == Y_INTERACTOR][['Structural_haplotype', 'ddG_REU']]
    m = x.merge(y, on='Structural_haplotype', suffixes=('_x', '_y'))
    m['Label'] = [s if isinstance(s, str) and s.strip() else str(h).split(',')[0].strip()
                  for s, h in zip(m['sHT'], m['Haplotypes'])]
    return m


def main():
    m = load()
    cross = (m['ddG_REU_x'].abs() > THRESHOLD) & (m['ddG_REU_y'].abs() > THRESHOLD)
    rest, hit = m[~cross], m[cross]
    strong = hit[hit['ddG_REU_y'].abs() > STRONG]
    plain = hit[hit['ddG_REU_y'].abs() <= STRONG]

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for t in XTICKS:
        ax.axvline(t, color=C_GRID, lw=0.6, zorder=0)
    for t in YTICKS:
        ax.axhline(t, color=C_GRID, lw=0.6, zorder=0)
    ax.axvline(0, color=C_ZERO, lw=0.8, ls=(0, (5, 4)), zorder=1)
    ax.axhline(0, color=C_ZERO, lw=0.8, ls=(0, (5, 4)), zorder=1)
    for t in (-THRESHOLD, THRESHOLD):
        ax.axvline(t, color=C_GRID, lw=0.7, ls=(0, (1, 2.5)), zorder=1)
        ax.axhline(t, color=C_GRID, lw=0.7, ls=(0, (1, 2.5)), zorder=1)

    ax.scatter(rest['ddG_REU_x'], rest['ddG_REU_y'], s=DOT_CLOUD, c=C_CLOUD,
               linewidths=0, alpha=ALPHA_CLOUD, zorder=2)
    ax.scatter(plain['ddG_REU_x'], plain['ddG_REU_y'], s=DOT_MARK, c=C_BLACK,
               linewidths=EDGE_MARK, zorder=3)
    for _, r in strong.iterrows():                       # one call per set, so each
        ax.scatter([r['ddG_REU_x']], [r['ddG_REU_y']],   # is its own vector group
                   s=DOT_MARK, c=COLOURS.get(r['Label'], C_BLACK),
                   linewidths=EDGE_MARK, zorder=4)

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks(XTICKS)
    ax.set_yticks(YTICKS)
    ax.set_xlabel(XLABEL, fontsize=FS_LAB)
    ax.set_ylabel(YLABEL, fontsize=FS_LAB)
    ax.tick_params(labelsize=FS_TICK, width=0.6, length=2.5, colors=C_AXIS)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(C_AXIS)

    os.makedirs(os.path.join(HERE, OUTDIR), exist_ok=True)
    fig.tight_layout()
    for f in FORMATS:
        fig.savefig(os.path.join(HERE, OUTDIR, 'figure2E.%s' % f), format=f,
                    dpi=PNG_DPI if f == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

    if WRITE_POSITIONS:
        out = hit[['Structural_haplotype', 'Label', 'ddG_REU_x', 'ddG_REU_y']].copy()
        out.columns = ['Structural_haplotype', 'Label', 'Receptor_ddG_REU', 'Pocket_factor_ddG_REU']
        out['Marker'] = np.where(out['Pocket_factor_ddG_REU'].abs() > STRONG,
                                 'coloured', 'black')
        out.sort_values('Receptor_ddG_REU').round(3) \
           .to_csv(os.path.join(HERE, OUTDIR, 'marked_sets_positions.csv'), index=False)
    print('written to %s' % OUTDIR)


if __name__ == '__main__':
    main()
