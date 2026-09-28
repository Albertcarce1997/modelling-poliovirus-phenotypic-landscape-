"""
Figure 1, panels F-Q: predicted changes in VP1 binding energy.

One panel per dataset and interactor.  Each VP1 structural haplotype is one dot; sets
are sorted by ddG up the vertical axis, so that axis is rank and the horizontal
axis is ddG.  Dashed lines mark the +/- 1 REU criterion.  All panels share one
symmetric ddG axis.  Reference sets in the Barcelona panels are ringed and named.

Run:      python figure1_panelsF-Q.py
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
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# ###########################################################################
# CONFIG
# ###########################################################################

OUTDIR  = 'figures'
FORMATS = ('pdf', 'svg', 'png')     # png is a preview only
PNG_DPI = 200

FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']

THRESHOLD = 1.0                     # REU; the criterion drawn as dashed lines

C_DESTAB = '#FF0000'                # ddG > +THRESHOLD
C_STAB   = '#00796B'                # ddG < -THRESHOLD
C_WITHIN = '#A6A6A6'                # within the criterion
C_AXIS   = '#333333'
C_GRID   = '#BFBFBF'

DOT_AUTO = True                     # scale the dot area to n
DOT      = 4.0                      # pt^2, used when DOT_AUTO is False
DOT_MIN, DOT_MAX = 1.6, 11.0
ALPHA    = 0.95

RANKPAD_FRAC = 0.025                # blank margin before the first and after the
                                    # last set, as a fraction of n

SCALE = 'global'                    # 'global' | 'interactor' | 'panel'
XLIM      = (-15.0, 15.0)           # symmetric, so zero sits at the centre
CLIP_MARK = True                    # sets outside XLIM are drawn at the edge
CLIP_SIZE = 9.0                     # area of that edge marker, pt^2
XTICK_N = 5                         # target ticks; the axis always ends on one

PANEL_W, PANEL_H = 2.05, 3.3        # inches per panel
FS_TICK, FS_LAB, FS_TITLE, FS_ANN = 7.0, 8.0, 8.5, 6.5

ANNOTATE_COUNTS = True
MARK_SETS    = True
LABEL_TEXT   = True
LABEL_HALO   = True                 # white halo so labels read over the dots
LABELS       = ['sHT-1', 'sHT-2', 'sHT-3', 'HT-21', 'NIE-ZAS']
LABEL_COLOUR = '#111111'
LABEL_PAD    = 5.0                  # marker-to-label clearance, in points
MARK_SIZE    = 14.0                 # open-circle area, pt^2
NOT_IN_N     = ['NIE-ZAS']          

DATASETS = ['Barcelona sewage (haplotype-resolved)',
            'Nigeria cVDPV2 outbreak (consensus)',
            'iVDPV2 (consensus)']
SHORT = {'Barcelona sewage (haplotype-resolved)': 'Barcelona',
         'Nigeria cVDPV2 outbreak (consensus)':   'Nigeria',
         'iVDPV2 (consensus)':                    'iVDPV2'}
LABEL_DATASETS = ['Barcelona sewage (haplotype-resolved)']

INTERACTORS = ['CD155 receptor', 'Pocket factor (palmitate)', 'mAb 9H2', 'mAb 10D2']
ITITLE = {'CD155 receptor': 'CD155 receptor', 'Pocket factor (palmitate)': 'Pocket factor',
          'mAb 9H2': 'mAb 9H2', 'mAb 10D2': 'mAb 10D2'}
ISLUG  = {'CD155 receptor': 'receptor', 'Pocket factor (palmitate)': 'pocket-factor',
          'mAb 9H2': '9H2', 'mAb 10D2': '10D2'}

COMBINED_PER_DATASET = True         # one four-panel figure per dataset
INDIVIDUAL_PANELS    = True         # and every panel on its own

# ###########################################################################
# end of CONFIG
# ###########################################################################

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('VDPV2_DATA') or os.path.join(
    HERE, 'VDPV2_VP1_binding_energy_raw_data.xlsx')
OUT = os.path.join(HERE, OUTDIR)

matplotlib.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': FONT,
    'pdf.fonttype': 42,             # text stays text in the PDF
    'svg.fonttype': 'none',         # and in the SVG
    'axes.linewidth': 0.6,
})
if matplotlib.font_manager.FontProperties(family=FONT).get_name() != FONT[0]:
    print('Note: %s was not found; the next available face in FONT is used.' % FONT[0])


def load():
    """One row per structural haplotype and interactor, with a display name."""
    d = pd.read_excel(DATA, 'Data')
    d = d[d['Structural_haplotype'] != 'Sabin 2 (reference)'].copy()

    def name(r):
        if isinstance(r['sHT'], str) and r['sHT'].strip():
            return r['sHT'].strip()
        h = str(r['Haplotypes'])
        return h.split(',')[0].strip() if h and h != 'nan' else str(r['Structural_haplotype'])

    d['Label'] = d.apply(name, axis=1)
    return d


def nice_limits(lo, hi, pad_it=True):
    """Round a range outward to whole ticks."""
    if pad_it:
        pad = 0.04 * (hi - lo)
        lo, hi = min(lo - pad, -THRESHOLD - pad), max(hi + pad, THRESHOLD + pad)
    ticks = MaxNLocator(nbins=XTICK_N, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = ticks[1] - ticks[0]
    lo2, hi2 = np.floor(lo / step) * step, np.ceil(hi / step) * step
    return lo2, hi2, np.arange(round(lo2 / step), round(hi2 / step) + 1) * step


def scale_limits(d):
    """ddG limits keyed by interactor; under 'global' every key is the same."""
    if XLIM is not None:
        lo, hi = XLIM
        steps = [st for st in (0.5, 1, 2, 2.5, 5, 10, 20) if
                 abs(round((hi - lo) / st) - (hi - lo) / st) < 1e-9]
        step = min(steps, key=lambda st: abs((hi - lo) / st - XTICK_N))
        lim = (lo, hi, np.arange(round(lo / step), round(hi / step) + 1) * step)
        return {i: lim for i in INTERACTORS}
    if SCALE == 'global':
        lim = nice_limits(d['ddG_REU'].min(), d['ddG_REU'].max())
        return {i: lim for i in INTERACTORS}
    if SCALE == 'interactor':
        return {i: nice_limits(g['ddG_REU'].min(), g['ddG_REU'].max())
                for i, g in d.groupby('Interactor')}
    return {}                                   # 'panel': fitted in panel()


def place_labels(ax, hit):
    """Name each marked set beside its circle, at the first free position
    around the marker: right, left, above, below, at increasing offsets."""
    bb = ax.get_window_extent()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x_per_pt = (x1 - x0) / max(bb.width / ax.figure.dpi * 72.0, 1.0)
    y_per_pt = (y1 - y0) / max(bb.height / ax.figure.dpi * 72.0, 1.0)
    th = 1.30 * FS_ANN * y_per_pt
    rx, ry = LABEL_PAD * x_per_pt, LABEL_PAD * y_per_pt

    mr = (MARK_SIZE / np.pi) ** 0.5 + 0.9          # marker radius, points
    boxes = [(r['ddG_REU'] - mr * x_per_pt, r['ddG_REU'] + mr * x_per_pt,
              r['rank'] - mr * y_per_pt, r['rank'] + mr * y_per_pt)
             for _, r in hit.iterrows()]
    out = []
    for _, r in hit.iterrows():
        w = len(str(r['Label'])) * FS_ANN * 0.58 * x_per_pt
        xv, yr = float(r['ddG_REU']), float(r['rank'])
        best = None
        for k in range(6):
            grow = 1.0 + 0.9 * k
            for pos in ('right', 'left', 'above', 'below'):
                if pos == 'right':
                    cx, cy, ha, va = xv + rx * grow, yr, 'left', 'center'
                    L, R, B, T = cx, cx + w, cy - th / 2, cy + th / 2
                elif pos == 'left':
                    cx, cy, ha, va = xv - rx * grow, yr, 'right', 'center'
                    L, R, B, T = cx - w, cx, cy - th / 2, cy + th / 2
                elif pos == 'above':
                    cx, cy, ha, va = xv, yr + ry * grow, 'center', 'bottom'
                    L, R, B, T = cx - w / 2, cx + w / 2, cy, cy + th
                else:
                    cx, cy, ha, va = xv, yr - ry * grow, 'center', 'top'
                    L, R, B, T = cx - w / 2, cx + w / 2, cy - th, cy
                if L < x0 or R > x1 or B < y0 or T > y1:
                    continue
                if any(not (R < l or L > rr or T < b or B > t) for l, rr, b, t in boxes):
                    continue
                best = (cx, cy, ha, va, L, R, B, T)
                break
            if best:
                break
        if best is None:
            cx, cy, ha, va = xv + rx, yr, 'left', 'center'
            best = (cx, cy, ha, va, cx, cx + w, cy - th / 2, cy + th / 2)
        cx, cy, ha, va, L, R, B, T = best
        boxes.append((L, R, B, T))
        out.append((cx, cy, ha, va, str(r['Label'])))

    for cx, cy, ha, va, txt in out:
        t = ax.annotate(txt, (cx, cy), ha=ha, va=va, fontsize=FS_ANN,
                        color=LABEL_COLOUR, zorder=7, clip_on=False)
        if LABEL_HALO:
            t.set_path_effects([pe.withStroke(linewidth=1.7, foreground='white'),
                                pe.Normal()])


def panel(ax, sub, dataset, xlim=None, title=None, xlabel=False, ylabel=False,
          report=None):
    s = sub.sort_values('ddG_REU', kind='mergesort').reset_index(drop=True)
    rank = np.arange(1, len(s) + 1)
    ddg = s['ddG_REU'].values
    col = np.where(ddg > THRESHOLD, C_DESTAB,
                   np.where(ddg < -THRESHOLD, C_STAB, C_WITHIN))

    pad = max(RANKPAD_FRAC * len(s), 0.6)
    ax.set_ylim(1 - pad, len(s) + pad)
    lo, hi, ticks = xlim if xlim else nice_limits(ddg.min(), ddg.max())
    ax.set_xlim(lo, hi)

    ax.axvline(0, color=C_AXIS, lw=0.6, zorder=2)
    for t in (THRESHOLD, -THRESHOLD):
        ax.axvline(t, color=C_GRID, lw=0.7, ls=(0, (3, 2)), zorder=1)

    sz = DOT
    if DOT_AUTO:
        hpt = ax.get_window_extent().height / ax.figure.dpi * 72.0
        sz = float(np.clip((0.75 * hpt / max(len(s), 1)) ** 2 * 3.0, DOT_MIN, DOT_MAX))
    inside = (ddg >= lo) & (ddg <= hi)
    ax.scatter(ddg[inside], rank[inside], s=sz, c=col[inside], linewidths=0,
               alpha=ALPHA, zorder=3)
    if CLIP_MARK:
        for beyond, at, mark in ((ddg > hi, hi, '>'), (ddg < lo, lo, '<')):
            if beyond.any():
                ax.scatter(np.full(beyond.sum(), at), rank[beyond], s=CLIP_SIZE,
                           c=col[beyond], marker=mark, linewidths=0,
                           clip_on=False, zorder=4)

    if MARK_SETS and dataset in LABEL_DATASETS:
        hit = s[s['Label'].isin(LABELS)].copy()
        hit['rank'] = hit.index + 1
        if len(hit):
            if LABEL_TEXT:
                place_labels(ax, hit)
            fill = np.where(hit['ddG_REU'] > THRESHOLD, C_DESTAB,
                            np.where(hit['ddG_REU'] < -THRESHOLD, C_STAB, C_WITHIN))
            ax.scatter(hit['ddG_REU'], hit['rank'], s=MARK_SIZE, c=fill,
                       edgecolors=LABEL_COLOUR, linewidths=0.7, zorder=8)
            if report is not None:
                for _, r in hit.sort_values('rank').iterrows():
                    report.append((r['Label'], int(r['rank']), round(float(r['ddG_REU']), 3)))

    ax.set_xticks(ticks)

    if ANNOTATE_COUNTS:
        dk = ddg[(~s['Label'].isin(NOT_IN_N)).values]
        rows = [('%d stabilising' % int((dk < -THRESHOLD).sum()), C_STAB),
                ('%d unchanged' % int((np.abs(dk) <= THRESHOLD).sum()), C_WITHIN),
                ('%d destabilising' % int((dk > THRESHOLD).sum()), C_DESTAB)]
        for k, (txt, colour) in enumerate(rows):
            ax.text(0.03, 0.985 - 0.038 * k, txt, transform=ax.transAxes,
                    va='top', ha='left', fontsize=FS_ANN, color=colour)

    n_sets = int((~s['Label'].isin(NOT_IN_N)).sum())
    if xlabel:
        ax.set_xlabel('ΔΔG (REU)', fontsize=FS_LAB)
    if ylabel:
        ax.set_ylabel('Structural haplotypes  (n = %d)' % n_sets,
                      fontsize=FS_LAB)
    if title:
        ax.set_title(title, fontsize=FS_TITLE, pad=4)
    ax.tick_params(labelsize=FS_TICK, width=0.6, length=2.5, colors=C_AXIS)
    ax.set_yticks([])
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color(C_AXIS)


def legend_handles(clipped=False):
    kw = dict(lw=0, marker='o', markersize=4)
    h = [Line2D([], [], color=C_STAB,
                label='\u0394\u0394G < \u2212%g REU (stabilising)' % THRESHOLD, **kw),
         Line2D([], [], color=C_WITHIN,
                label='within \u00b1%g REU' % THRESHOLD, **kw),
         Line2D([], [], color=C_DESTAB,
                label='\u0394\u0394G > +%g REU (destabilising)' % THRESHOLD, **kw)]
    if clipped:
        h.append(Line2D([], [], color=C_DESTAB, lw=0, marker='>', markersize=4,
                        label='beyond the axis'))
    return h



def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    for f in FORMATS:
        fig.savefig(os.path.join(OUT, '%s.%s' % (stem, f)), format=f,
                    dpi=PNG_DPI if f == 'png' else None,
                    bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


def main():
    d = load()
    lims = scale_limits(d)
    rep = []
    for ds in DATASETS:
        sub_ds = d[d['Dataset'] == ds]
        if COMBINED_PER_DATASET:
            fig, axes = plt.subplots(1, len(INTERACTORS),
                                     figsize=(PANEL_W * len(INTERACTORS), PANEL_H))
            for k, it in enumerate(INTERACTORS):
                r = []
                panel(axes[k], sub_ds[sub_ds['Interactor'] == it], ds,
                      xlim=lims.get(it), title=ITITLE[it], xlabel=True,
                      ylabel=(k == 0), report=r)
                rep += [(SHORT[ds], it) + t for t in r]
            clipped = bool(CLIP_MARK and
                           ((sub_ds['ddG_REU'] > lims[INTERACTORS[0]][1]).any() or
                            (sub_ds['ddG_REU'] < lims[INTERACTORS[0]][0]).any()))
            fig.legend(handles=legend_handles(clipped), loc='lower center', ncol=4,
                       frameon=False, fontsize=FS_ANN, bbox_to_anchor=(0.5, -0.065))
            fig.suptitle(SHORT[ds], fontsize=FS_TITLE + 1, y=1.02)
            fig.tight_layout()
            save(fig, 'figure1_%s' % SHORT[ds])
        if INDIVIDUAL_PANELS:
            for it in INTERACTORS:
                fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
                panel(ax, sub_ds[sub_ds['Interactor'] == it], ds,
                      xlim=lims.get(it), title=ITITLE[it], xlabel=True, ylabel=True)
                fig.tight_layout()
                save(fig, 'figure1_%s_%s' % (SHORT[ds], ISLUG[it]))
    if rep:
        pd.DataFrame(rep, columns=['Dataset', 'Interactor', 'Set', 'Rank', 'ddG_REU']) \
          .to_csv(os.path.join(OUT, 'marked_sets_positions.csv'), index=False)
    print('written to %s' % OUT)


if __name__ == '__main__':
    main()
