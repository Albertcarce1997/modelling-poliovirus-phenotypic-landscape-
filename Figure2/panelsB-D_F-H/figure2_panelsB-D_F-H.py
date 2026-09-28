"""
Figure 2, panels B-D and F-H: interface differences against Sabin 2.

One panel per interactor and measure.  Each haplotype is a dot at its mean
difference from Sabin 2, with a bar from mean - SD to mean + SD, SD being the
standard deviation across the models of its ensemble.  Colour marks the side of
zero.  Sabin 2 is the top row of every panel, grey and at zero by definition, so
its own spread between models can be read against the haplotypes below it.  All
panels of one measure share one symmetric axis.

Run:      python figure2_panelsB-D_F-H.py
Data:     Figure2_interface_descriptors.xlsx  (or $VDPV2_FIG2_DATA)
Output:   figures/, as PDF, SVG and PNG with text kept as text
Requires: python >= 3.9, matplotlib >= 3.6, numpy, pandas, openpyxl
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ###########################################################################
# CONFIG
# ###########################################################################

OUTDIR  = 'figures'
FORMATS = ('pdf', 'svg', 'png')     # png is a preview only
PNG_DPI = 200

FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']

REFERENCE = 'Sabin 2'               # the group every other one is measured from

C_UP     = '#D55E00'                # above the reference
C_DOWN   = '#0072B2'                # below the reference
C_FLAT   = '#A6A6A6'                # a haplotype that does not differ
C_REF    = '#4D4D4D'                # the reference row
C_BAND   = '#D9D9D9'                # +/- 1 SD of the reference, when REF_BAND
C_AXIS   = '#333333'
C_GRID   = '#BFBFBF'
COLOUR_BY_SIGN = True               # False paints every marker C_AXIS

DOT      = 22.0                     # dot area, pt^2
SD_BARS  = True                     # mean +/- 1 SD of the haplotype, across models
SD_LW    = 1.2
SD_CAP   = 2.0                      # cap half-height, in points
REF_ROW  = True                     # the reference as the top row of each panel
REF_GAP  = 0.4                      # its clearance from the haplotypes, in rows
REF_LABEL_VALUE = False             # its value is zero by definition
REF_BAND = False                    # +/- 1 SD of the reference, shaded round zero

VALUE_LABELS = True
VALUE_FMT = {'BSA': '%+.1f', 'Contacting residues': '%+.2f', 'H-bonds': '%+.2f'}
VALUE_PAD = 3.0                     # dot-to-label clearance, in points

ORDER = 'data'                      # 'data' keeps the order of the violin plots;
                                    # 'value' sorts each panel by its difference
                                    # and the panels then no longer align
ROWPAD_FRAC = 0.6                   # blank margin above the first and below the
                                    # last haplotype, in row units

SCALE = 'measure'                   # 'measure' shares one axis across every panel
                                    # of a measure; 'panel' gives each panel its own

# one symmetric axis per measure; None computes it from the data
XLIM = {'BSA': (-150.0, 150.0),
        'Contacting residues': (-6.0, 6.0),
        'H-bonds': (-4.0, 4.0)}
XTICK_STEP = {'BSA': 50.0, 'Contacting residues': 2.0, 'H-bonds': 2.0}
XTICK_N = 5                         # target ticks when XLIM is None

PANEL_W, PANEL_H = 2.15, 2.75       # inches per panel
FS_TICK, FS_LAB, FS_TITLE, FS_ANN = 7.0, 8.0, 8.5, 6.5

PARTS = ['2A-D', '2E-H']
MEASURES = ['BSA', 'Contacting residues', 'H-bonds']
INTERACTORS = ['CD155 receptor', 'Pocket factor (palmitate)', 'mAb 9H2', 'mAb 10D2']

AXLABEL = {'BSA': 'Difference in BSA (Å²)',
           'Contacting residues': 'Difference in contacting residues',
           'H-bonds': 'Difference in hydrogen bonds'}
ITITLE = {'CD155 receptor': 'CD155 receptor', 'Pocket factor (palmitate)': 'Pocket factor',
          'mAb 9H2': 'mAb 9H2', 'mAb 10D2': 'mAb 10D2'}
ISLUG  = {'CD155 receptor': 'receptor', 'Pocket factor (palmitate)': 'pocket-factor',
          'mAb 9H2': '9H2', 'mAb 10D2': '10D2'}
MSLUG  = {'BSA': 'BSA', 'Contacting residues': 'residues', 'H-bonds': 'hbonds'}
PSLUG  = {'2A-D': 'panelsB-D', '2E-H': 'panelsF-H'}
PTITLE = {'2A-D': 'Panels 2B\u20132D', '2E-H': 'Panels 2F\u20132H'}

COMBINED_PER_PART = True            # one grid of every panel of the part
ROW_PER_MEASURE   = True            # one row of interactors per measure
INDIVIDUAL_PANELS = True            # and every panel on its own

# ###########################################################################
# end of CONFIG
# ###########################################################################

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': FONT,
                     'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 0.6})


HERE = os.path.dirname(os.path.abspath(__file__))


def _find_data():
    return os.environ.get('VDPV2_FIG2_DATA') or os.path.join(HERE, 'Figure2_interface_descriptors.xlsx')


def load():
    d = pd.read_excel(_find_data(), sheet_name='Data')
    for c in ('Part', 'Panel', 'Measure', 'Interactor', 'Group'):
        d[c] = d[c].astype(str).str.strip()
    d['Value'] = pd.to_numeric(d['Value'], errors='coerce')
    return d.dropna(subset=['Value'])


def summarise(d):
    """Group means and SDs, and each group's difference from the reference."""
    g = (d.groupby(['Part', 'Panel', 'Measure', 'Interactor', 'Group'])['Value']
          .agg(n='size', mean='mean',
               sd=lambda v: v.std(ddof=1) if len(v) > 1 else 0.0)
          .reset_index())
    ref = (g[g['Group'] == REFERENCE]
           .set_index(['Part', 'Measure', 'Interactor'])[['mean', 'sd', 'n']])
    key = list(zip(g['Part'], g['Measure'], g['Interactor']))
    g['ref_mean'] = [ref['mean'].get(k, np.nan) for k in key]
    g['ref_sd'] = [ref['sd'].get(k, np.nan) for k in key]
    g['ref_n'] = [ref['n'].get(k, np.nan) for k in key]
    g['diff'] = g['mean'] - g['ref_mean']
    return g[g['Group'] != REFERENCE].copy()


def group_order(d, part):
    seen = [g for g in dict.fromkeys(d[d['Part'] == part]['Group']) if g != REFERENCE]
    return seen


def scale_limits(s):
    """One symmetric limit and tick set per measure, or per panel."""
    keys = ([(m,) for m in MEASURES] if SCALE == 'measure'
            else [(m, i) for m in MEASURES for i in INTERACTORS])
    lims = {}
    for k in keys:
        m = k[0]
        sub = s[s['Measure'] == m]
        if len(k) > 1:
            sub = sub[sub['Interactor'] == k[1]]
        if sub.empty:
            continue
        if XLIM and m in XLIM and XLIM[m] is not None:
            lo, hi = XLIM[m]
            step = (XTICK_STEP or {}).get(m)
            if not step:
                step = (hi - lo) / max(2, XTICK_N - 1)
        else:
            span = float(np.nanmax(np.abs(np.concatenate(
                [(sub['diff'] - sub['sd']).values, (sub['diff'] + sub['sd']).values]))))
            step = _nice_step(2 * span / max(2, XTICK_N - 1))
            hi = step * np.ceil(span / step)
            lo = -hi
        ticks = np.arange(round(lo / step), round(hi / step) + 1) * step
        lims[k] = (lo, hi, ticks)
    return lims


def limit_for(lims, measure, interactor):
    return lims[(measure,)] if SCALE == 'measure' else lims[(measure, interactor)]


def _nice_step(raw):
    if raw <= 0:
        return 1.0
    e = 10.0 ** np.floor(np.log10(raw))
    for f in (1, 2, 2.5, 5, 10):
        if raw <= f * e:
            return f * e
    return 10 * e


def panel(ax, s, order, measure, lim, title=None, xlabel=True, ylabel=True):
    s = s.set_index('Group').reindex(order)
    if ORDER == 'value':
        s = s.sort_values('diff')
    ref_sd = float(s['ref_sd'].iloc[0]) if np.isfinite(s['ref_sd'].iloc[0]) else 0.0
    labels = list(s.index)
    dv = s['diff'].values
    sd = s['sd'].values
    is_ref = np.zeros(len(s), bool)
    if REF_ROW:                                    # the reference joins as a row
        labels = [REFERENCE] + labels
        dv = np.concatenate([[0.0], dv])
        sd = np.concatenate([[ref_sd], sd])
        is_ref = np.concatenate([[True], is_ref])
    y = np.arange(len(labels))[::-1].astype(float)  # first row at the top
    if REF_ROW:
        y[0] += REF_GAP
    lo, hi, ticks = lim

    if REF_BAND and np.isfinite(s['ref_sd'].iloc[0]) and s['ref_sd'].iloc[0] > 0:
        r = float(s['ref_sd'].iloc[0])
        ax.axvspan(-r, r, color=C_BAND, lw=0, zorder=0)
    ax.axvline(0.0, color=C_AXIS, lw=0.8, zorder=2)
    ax.set_xticks(ticks)
    ax.xaxis.grid(True, color=C_GRID, lw=0.4, alpha=0.6, zorder=1)
    ax.set_axisbelow(True)

    if COLOUR_BY_SIGN:
        col = np.where(dv > 0, C_UP, np.where(dv < 0, C_DOWN, C_FLAT))
    else:
        col = np.full(len(dv), C_AXIS)
    col = np.where(is_ref, C_REF, col)

    if SD_BARS:
        for yi, v, e, c in zip(y, dv, sd, col):
            if not np.isfinite(e) or e <= 0:
                continue
            ax.plot([v - e, v + e], [yi, yi], color=c, lw=SD_LW, zorder=4)
            for x in (v - e, v + e):
                ax.plot([x, x], [yi, yi], color=c, marker='|', markersize=2 * SD_CAP,
                        mew=SD_LW, zorder=4)
    ax.scatter(dv, y, s=DOT, c=col, edgecolors='white', linewidths=0.5, zorder=5)

    if VALUE_LABELS:
        fmt = VALUE_FMT.get(measure, '%+.2f')
        for yi, v, e, r in zip(y, dv, sd, is_ref):
            if r and not REF_LABEL_VALUE:
                continue
            end = v + (e if v >= 0 else -e) if (SD_BARS and np.isfinite(e)) else v
            right = v >= 0
            ax.annotate(fmt % v, (end, yi), textcoords='offset points',
                        xytext=(VALUE_PAD if right else -VALUE_PAD, 0),
                        ha='left' if right else 'right', va='center',
                        fontsize=FS_ANN, color=C_AXIS, annotation_clip=False, zorder=6)

    ax.set_xlim(lo, hi)
    ax.set_ylim(y.min() - ROWPAD_FRAC, y.max() + ROWPAD_FRAC)
    ax.set_yticks(y)
    ax.set_yticklabels(labels if ylabel else [])
    if ylabel and REF_ROW:
        ax.get_yticklabels()[0].set_color(C_REF)
    if xlabel:
        ax.set_xlabel(AXLABEL.get(measure, 'Difference'), fontsize=FS_LAB)
    if title:
        ax.set_title(title, fontsize=FS_TITLE, pad=4)
    ax.tick_params(labelsize=FS_TICK, width=0.6, length=2.5, colors=C_AXIS)
    ax.tick_params(axis='y', length=0)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color(C_AXIS)


def legend_handles():
    h = []
    if COLOUR_BY_SIGN:
        h += [Line2D([], [], color=C_UP, lw=0, marker='o', markersize=4,
                     label='above %s' % REFERENCE),
              Line2D([], [], color=C_DOWN, lw=0, marker='o', markersize=4,
                     label='below %s' % REFERENCE),
              Line2D([], [], color=C_FLAT, lw=0, marker='o', markersize=4,
                     label='no difference')]
    else:
        h += [Line2D([], [], color=C_AXIS, lw=0, marker='o', markersize=4,
                     label='difference from %s' % REFERENCE)]
    if REF_ROW:
        h.append(Line2D([], [], color=C_REF, lw=0, marker='o', markersize=4,
                        label='%s (reference)' % REFERENCE))
    h.append(Line2D([], [], color=C_AXIS, lw=SD_LW, marker='|', markersize=2 * SD_CAP,
                    label='mean ± 1 SD across models'))
    if REF_BAND:
        h.append(Patch(facecolor=C_BAND, edgecolor='none',
                       label='± 1 SD of %s' % REFERENCE))
    return h


def save(fig, stem, formats=FORMATS):
    os.makedirs(os.path.join(HERE, OUTDIR), exist_ok=True)
    for f in formats:
        fig.savefig(os.path.join(HERE, OUTDIR, '%s.%s' % (stem, f)), format=f,
                    dpi=PNG_DPI if f == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


def main():
    d = load()
    s = summarise(d)
    lims = scale_limits(s)

    for part in PARTS:
        sp = s[s['Part'] == part]
        if sp.empty:
            continue
        order = group_order(d, part)
        inters = [i for i in INTERACTORS if i in set(sp['Interactor'])]
        meas = [m for m in MEASURES if m in set(sp['Measure'])]

        if COMBINED_PER_PART:
            fig, axes = plt.subplots(len(meas), len(inters), squeeze=False,
                                     figsize=(PANEL_W * len(inters) + 0.5,
                                              PANEL_H * len(meas)))
            for r, m in enumerate(meas):
                for c, it in enumerate(inters):
                    panel(axes[r][c],
                          sp[(sp['Measure'] == m) & (sp['Interactor'] == it)],
                          order, m, limit_for(lims, m, it),
                          title=ITITLE[it] if r == 0 else None,
                          xlabel=True, ylabel=(c == 0))
            fig.legend(handles=legend_handles(), loc='lower center',
                       ncol=5 if len(inters) >= 3 else 3,
                       frameon=False, fontsize=FS_ANN,
                       bbox_to_anchor=(0.5, -0.02 - 0.01 * len(meas)))
            fig.suptitle(PTITLE[part], fontsize=FS_TITLE + 1, y=1.0)
            fig.tight_layout()
            save(fig, 'figure2_%s_all' % PSLUG[part])

        if ROW_PER_MEASURE:
            for m in meas:
                fig, axes = plt.subplots(1, len(inters), squeeze=False,
                                         figsize=(PANEL_W * len(inters) + 0.5, PANEL_H))
                for c, it in enumerate(inters):
                    panel(axes[0][c],
                          sp[(sp['Measure'] == m) & (sp['Interactor'] == it)],
                          order, m, limit_for(lims, m, it), title=ITITLE[it],
                          xlabel=True, ylabel=(c == 0))
                fig.legend(handles=legend_handles(), loc='lower center',
                           ncol=5 if len(inters) >= 3 else 3, frameon=False,
                           fontsize=FS_ANN, bbox_to_anchor=(0.5, -0.08))
                fig.tight_layout()
                save(fig, 'figure2_%s_%s' % (PSLUG[part], MSLUG[m]))

        if INDIVIDUAL_PANELS:
            for m in meas:
                for it in inters:
                    fig, ax = plt.subplots(figsize=(PANEL_W + 0.5, PANEL_H))
                    panel(ax, sp[(sp['Measure'] == m) & (sp['Interactor'] == it)],
                          order, m, limit_for(lims, m, it), title=ITITLE[it],
                          xlabel=True, ylabel=True)
                    fig.tight_layout()
                    save(fig, 'figure2_%s_%s_%s' % (PSLUG[part], MSLUG[m], ISLUG[it]),
                         formats=tuple(f for f in FORMATS if f != 'png'))

    os.makedirs(os.path.join(HERE, OUTDIR), exist_ok=True)
    out = s[['Part', 'Panel', 'Measure', 'Interactor', 'Group', 'n', 'mean', 'sd',
             'ref_n', 'ref_mean', 'ref_sd', 'diff']].copy()
    out.columns = ['Part', 'Panel', 'Measure', 'Interactor', 'Group', 'n_group',
                   'Mean_group', 'SD_group', 'n_reference', 'Mean_reference',
                   'SD_reference', 'Difference']
    out.round(4).to_csv(os.path.join(HERE, OUTDIR, 'plotted_differences.csv'), index=False)
    print('written to %s' % OUTDIR)


if __name__ == '__main__':
    main()
