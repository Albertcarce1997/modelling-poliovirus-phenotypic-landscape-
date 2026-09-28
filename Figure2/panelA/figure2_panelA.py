"""
Figure 2, panel A - the two 3-D VDPV2 VP1 maps, as editable vector graphics.

Markers are small and semi-transparent (MARKER_ALPHA), so that overlapping
points in the dense part of the cloud remain distinguishable.

Plain black markers, grid and axes only: no colours, ellipses, leader lines or
haplotype labels, so that annotation can be done in Illustrator on the vector
output.  The highlighted haplotypes are drawn as a SEPARATE scatter on top of
the cloud, so in Illustrator they arrive as their own group and are easy to
select and recolour one by one.

Run:  python figure2_panelA.py            (writes ./figures)
Data: ../data/Figure2A_ddG_by_mutation_set.xlsx, ../data/Figure2A_haplotypes_per_set.xlsx

Requires: matplotlib >= 3.8, numpy, pandas, openpyxl   (pypdf only for COMBINE)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
import pandas as pd

# ###########################################################################
# CONFIG
# ###########################################################################

DATA_FILE = 'VDPV2_VP1_binding_energy_raw_data.xlsx'
DATASET = 'Barcelona sewage (haplotype-resolved)'

OUTDIR = 'figures'
FORMATS = ('pdf', 'svg', 'png')   # png is only a preview, for iterating
PNG_DPI = 200
COMBINE = True               # also write both panels on one page (needs pypdf)

# ---------------------------------------------------------------- appearance
FONT = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']
MARKER_SIZE = 12              # area in pt^2 of the mutation-set markers
HIGHLIGHT_SIZE = 15          # area in pt^2 of the highlighted haplotypes
MARKER_COLOUR = '#111111'
MARKER_ALPHA = 0.65          # marker transparency (0 = clear, 1 = opaque)
GRID_COLOUR = '#CCCCCC'      # gridlines, pane edges and axis lines
GRID_LW = 0.7
MINUS = '−'             # true minus sign in tick labels

# Structural haplotypes to leave out at plotting time. 
EXCLUDE_SETS = []

# Haplotypes shown in the zoom map. 
COLOUR_HIGHLIGHTED = True
HIGHLIGHTED = {
    'HT-731': '#B03F1F',
    'HT-750': '#A6D948',
    'HT-325': '#C281C7',   
    'HT-626': '#49BB39',
    'HT-987': '#939149',
    'HT-257': '#413A16',
    'HT-422': '#2B8D3C',
    'HT-294': '#5DA2DD',
    'HT-346': '#237335',
    'HT-26':  '#396EBC',
    'HT-148': '#7B6F83',   
    'HT-69':  '#F3E7AF',
    'HT-518': '#EC078F',
}

# --------------------------------------------------------------------- panels
# view      : (elevation, azimuth, perspective focal length) changing these changes the orientation.
# box       : axes box aspect (x, y, z).
# lim       : axis limits.  x = 9H2 ddG, y = -(receptor ddG)  (so that the
#             receptor axis reads 3 -> -3), z = 10D2 ddG.
# trim      : shrink the far end of each axis until the box just encloses the
#             ticks and the data (+ pad).  The markers do NOT move and the
#             perspective is unchanged.
#             Set False to keep the cropped box.
# trim_pad  : slack left beyond the outermost tick / point, in data units.
# ticks     : tick values per axis (y in plotted units, i.e. -receptor).
# tick_edge : which box edge each axis' ticks sit on, as (side, side) for the
#             two OTHER axes: 'lo' = low limit, 'hi' = high limit.
#             x -> (y side, z side);  y -> (x side, z side);  z -> (x side, y side)
# tick_off  : tick label offset from the box edge, in points (x right, y up).
# label_off : axis title offset from the middle of the same edge, in points.
# scale_ref : (haplotype A, haplotype B, distance in pt) - sets the plotting
#             scale.  
# figsize   : canvas (width, height) in points, or None to size it so that
#             everything drawn fits exactly.
# centre_on : what to centre in the canvas when figsize is given -
#             'all' = everything drawn, 'data' = the markers only.
# margin    : white space kept around everything drawn, in points.
#
# Common tweaks
#   box does not close / something is cut off -> figsize=None, or raise margin
#   a tick label sits on the grid            -> change that axis' tick_off
#   an axis title sits on the grid           -> change that axis' label_off
#   ticks on the wrong edge of the box       -> change that axis' tick_edge
#   too much empty box below the cloud       -> figsize=(W, H) + centre_on='data'
#   whole panel too large / too small        -> scale both scale_ref numbers
#                                               and MARKER_SIZE by one factor

OVERVIEW = dict(
    name='Figure2A_overview',
    which='all',                       # all structural haplotypes
    view=(0.492, -31.156, 0.198),
    box=(1.0, 0.689, 0.832),
    lim=dict(x=(-10.324, 15.305), y=(-3.952, 8.043), z=(-10.324, 15.305)),
    trim=True,
    trim_pad=dict(x=0.45, y=0.30, z=0.45),
    ticks=dict(x=[-10, -5, 0, 5, 10],
               y=[-3, -2, -1, 0, 1, 2, 3],
               z=[-10, -5, 0, 5, 10]),
    tick_edge=dict(x=('lo', 'lo'), y=('hi', 'lo'), z=('lo', 'lo')),
    tick_off=dict(x=(-11.0, -1.8), y=(-13.0, 10.4), z=(-23.0, 0.0)),
    label=dict(x='9H2 mAb ΔΔG (REU)',
               y='Receptor ΔΔG (REU)',
               z='10D2 mAb ΔΔG (REU)'),
    label_off=dict(x=(-30.0, -14.0), y=(14.0, 3.0), z=(-46.0, 0.0)),
    tick_size=9.0,
    label_size=10.5,
    scale_ref=('HT-518', 'HT-422', 142.40),
    figsize=None,
    centre_on='all',
    margin=8.0,
)

ZOOM = dict(
    name='Figure2A_zoom',
    which='highlighted',               # only the 13 haplotypes
    view=(15.849, -52.783, 0.247),
    box=(1.0, 1.0, 1.0),
    lim=dict(x=(0.670, 14.141), y=(1.032, 3.198), z=(0.766, 7.055)),
    trim=True,
    trim_pad=dict(x=0.30, y=0.10, z=0.15),
    ticks=dict(x=[2, 4, 6, 8, 10, 12, 14],
               y=[1.5, 2.0, 2.5, 3.0],
               z=[1, 2, 3, 4, 5, 6, 7]),
    tick_edge=dict(x=('lo', 'lo'), y=('hi', 'lo'), z=('lo', 'lo')),
    tick_off=dict(x=(-9.0, -5.0), y=(9.0, -5.0), z=(-12.0, 0.0)),
    label=dict(x='9H2 mAb ΔΔG (REU)',
               y='Receptor ΔΔG (REU)',
               z='10D2 mAb ΔΔG (REU)'),
    label_off=dict(x=(-22.0, -18.0), y=(20.0, -20.0), z=(-30.0, 0.0)),
    tick_size=7.5,
    label_size=9.0,
    scale_ref=('HT-26', 'HT-422', 169.95),
    figsize=None,
    centre_on='all',
    margin=8.0,
)

PANELS = [OVERVIEW, ZOOM]

# ###########################################################################
# end of CONFIG
# ###########################################################################

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('VDPV2_DATA') or os.path.join(HERE, DATA_FILE)
OUT = os.path.join(HERE, OUTDIR)
os.makedirs(OUT, exist_ok=True)

matplotlib.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': FONT,
    'pdf.fonttype': 42,          # keep text as text in the PDF
    'svg.fonttype': 'none',      # ditto in the SVG
    'axes.unicode_minus': True,
})


def load():
    """The Barcelona structural haplotypes as (x, y, z) = (9H2, -receptor, 10D2), plus
       the row index of every haplotype name."""
    d = pd.read_excel(DATA, sheet_name='Data')
    d = d[(d['Dataset'] == DATASET) & (d['Structural_haplotype'].str.startswith('Set_'))]
    wide = d.pivot_table(index='Structural_haplotype', columns='Interactor', values='ddG_REU')
    wide = wide.reindex(sorted(wide.index, key=lambda s: int(s.split('_')[1])))
    P = np.column_stack([wide['mAb 9H2'].to_numpy(float),
                         -wide['CD155 receptor'].to_numpy(float),
                         wide['mAb 10D2'].to_numpy(float)])
    keep = np.array([m not in EXCLUDE_SETS for m in wide.index])

    names = (d[d['Interactor'] == 'CD155 receptor']
             .set_index('Structural_haplotype')['Haplotypes'].reindex(wide.index))
    idx = {}
    for i, entry in enumerate(names):
        if isinstance(entry, str):
            for h in entry.split(','):
                idx[h.strip()] = i
    return P, idx, keep


def fmt(v):
    return ('%g' % v).replace('-', MINUS)


def proj(ax, pts):
    """3-D points -> display coordinates"""
    M = ax.get_proj()
    return np.array([ax.transData.transform(proj3d.proj_transform(*p, M)[:2])
                     for p in pts])


def trim_box(cfg, P, sel):
    """Cut the unused far end off each axis, adjusting the box aspect by the
       same factor so that the projection of the data is untouched."""
    lim = {k: list(v) for k, v in cfg['lim'].items()}
    box = list(cfg['box'])
    if not cfg.get('trim'):
        return lim, box
    for i, k in enumerate('xyz'):
        lo, hi = lim[k]
        need = max(max(cfg['ticks'][k]), P[sel, i].max()) + cfg['trim_pad'][k]
        s = (need - lo) / (hi - lo)
        if s < 1.0:
            lim[k][1] = lo + s * (hi - lo)
            box[i] *= s
    return lim, box


def build(cfg, P, idx, keep):
    sel = np.flatnonzero(keep) if cfg['which'] == 'all' else \
        np.array([idx[h] for h in HIGHLIGHTED])
    hi_sel = [idx[h] for h in HIGHLIGHTED]
    lim, box = trim_box(cfg, P, sel)

    fig = plt.figure(figsize=(6, 6))          # resized once the layout is known
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(*lim['x']); ax.set_ylim(*lim['y']); ax.set_zlim(*lim['z'])
    ax.set_box_aspect(box)
    ax.set_proj_type('persp', focal_length=cfg['view'][2])
    ax.view_init(elev=cfg['view'][0], azim=cfg['view'][1])
    ax.set_xticks(cfg['ticks']['x'])
    ax.set_yticks(cfg['ticks']['y'])
    ax.set_zticks(cfg['ticks']['z'])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor(GRID_COLOUR); axis.pane.set_linewidth(GRID_LW)
        axis.line.set_color(GRID_COLOUR); axis.line.set_linewidth(GRID_LW)
        axis._axinfo['grid'].update(color=GRID_COLOUR, linewidth=GRID_LW)
        axis.set_tick_params(colors='none', labelsize=0.01, length=0)
    ax.set_xlabel(''); ax.set_ylabel(''); ax.set_zlabel('')
    ax.grid(True)

    rest = np.setdiff1d(sel, np.array(hi_sel))
    if len(rest):
        ax.scatter(P[rest, 0], P[rest, 1], P[rest, 2], s=MARKER_SIZE,
                   c=MARKER_COLOUR, alpha=MARKER_ALPHA, depthshade=False,
                   linewidths=0)
    keep = [i for i in hi_sel if i in set(sel.tolist())]
    if keep:                           # separate group, drawn on top
        cols = [HIGHLIGHTED[h] if COLOUR_HIGHLIGHTED else MARKER_COLOUR
                for h in HIGHLIGHTED if idx[h] in keep]
        order = [idx[h] for h in HIGHLIGHTED if idx[h] in keep]
        ax.scatter(P[order, 0], P[order, 1], P[order, 2], s=HIGHLIGHT_SIZE,
                   c=cols, alpha=MARKER_ALPHA, depthshade=False, linewidths=0)
    return fig, ax, lim, box, sel


def edge_points(cfg, lim, which, values):
    """3-D positions of `values` along `which`, on the configured box edge"""
    a, b = cfg['tick_edge'][which]
    others = {'x': ('y', 'z'), 'y': ('x', 'z'), 'z': ('x', 'y')}[which]
    p = lim[others[0]][0 if a == 'lo' else 1]
    q = lim[others[1]][0 if b == 'lo' else 1]
    return [{'x': (v, p, q), 'y': (p, v, q), 'z': (p, q, v)}[which] for v in values]


def layout(fig, ax, cfg, lim, P, idx, sel):
    """Scale the axes so the reference distance matches, then size and centre
       the canvas around everything that gets drawn."""
    a, b, target = cfg['scale_ref']
    ref = [P[idx[a]], P[idx[b]]]
    for _ in range(3):                 # linear, so this converges at once
        fig.canvas.draw()
        d = proj(ax, ref) * 72.0 / fig.dpi
        cur = np.hypot(*(d[1] - d[0]))
        if abs(cur - target) < 1e-6:
            break
        s = target / cur
        p = ax.get_position()
        ax.set_position([p.x0, p.y0, p.width * s, p.height * s])

    # everything that will be inked, in points
    fig.canvas.draw()
    pts = [(x, y, z) for x in lim['x'] for y in lim['y'] for z in lim['z']]
    D = list(proj(ax, pts) * 72.0 / fig.dpi)
    for w in 'xyz':
        E = proj(ax, edge_points(cfg, lim, w, cfg['ticks'][w])) * 72.0 / fig.dpi
        o = np.array(cfg['tick_off'][w])
        half = 0.62 * cfg['tick_size'] * np.array([2.2, 1.0])   # rough text box
        D += [e + o + half for e in E] + [e + o - half for e in E]
        M = proj(ax, edge_points(cfg, lim, w,
                                 [cfg['ticks'][w][0], cfg['ticks'][w][-1]]))
        M = M.mean(0) * 72.0 / fig.dpi + np.array(cfg['label_off'][w])
        r = 0.5 * cfg['label_size'] * len(cfg['label'][w]) * 0.55
        D += [M + r, M - r]
    D = np.array(D)
    m = cfg['margin']
    lo, hi = D.min(0) - m, D.max(0) + m
    if cfg['figsize']:
        W, H = cfg['figsize']
        if cfg.get('centre_on') == 'data':      # crop around the markers
            Q = proj(ax, P[sel]) * 72.0 / fig.dpi
            lo, hi = Q.min(0) - m, Q.max(0) + m
    else:
        W, H = hi - lo

    p = ax.get_position()
    fw, fh = fig.get_size_inches() * 72.0
    x0, y0 = p.x0 * fw, p.y0 * fh
    shift = (np.array([W, H]) - (hi - lo)) / 2.0 - lo
    fig.set_size_inches(W / 72.0, H / 72.0)
    ax.set_position([(x0 + shift[0]) / W, (y0 + shift[1]) / H,
                     p.width * fw / W, p.height * fh / H])
    fig.canvas.draw()
    d = proj(ax, ref) * 72.0 / fig.dpi
    return W, H, np.hypot(*(d[1] - d[0]))


def annotate(fig, ax, cfg, lim):
    W, H = fig.get_size_inches() * 72.0
    ov = fig.add_axes([0, 0, 1, 1], zorder=5)
    ov.set_xlim(0, W); ov.set_ylim(0, H); ov.axis('off'); ov.patch.set_alpha(0)
    for w in 'xyz':
        vals = cfg['ticks'][w]
        D = proj(ax, edge_points(cfg, lim, w, vals)) * 72.0 / fig.dpi
        o = cfg['tick_off'][w]
        for (x, y), v in zip(D, vals):
            ov.text(x + o[0], y + o[1], fmt(-v if w == 'y' else v),
                    ha='center', va='center', fontsize=cfg['tick_size'],
                    color='black')
        E = proj(ax, edge_points(cfg, lim, w, [vals[0], vals[-1]])) * 72.0 / fig.dpi
        mid = E.mean(0)
        dx, dy = E[1] - E[0]
        rot = np.degrees(np.arctan2(dy, dx))
        rot = rot - 180 if rot > 90 else (rot + 180 if rot < -90 else rot)
        lo_ = cfg['label_off'][w]
        ov.text(mid[0] + lo_[0], mid[1] + lo_[1], cfg['label'][w], rotation=rot,
                rotation_mode='anchor', ha='center', va='center',
                fontsize=cfg['label_size'], color='black')
    return ov


def save(fig, name):
    paths = []
    for ext in FORMATS:
        p = os.path.join(OUT, '%s.%s' % (name, ext))
        fig.savefig(p, transparent=(ext != 'png'),
                    dpi=PNG_DPI if ext == 'png' else None,
                    facecolor='white' if ext == 'png' else 'none')
        paths.append(p)
    return paths


def main():
    P, idx, keep = load()
    for cfg in PANELS:
        fig, ax, lim, box, sel = build(cfg, P, idx, keep)
        W, H, got = layout(fig, ax, cfg, lim, P, idx, sel)
        annotate(fig, ax, cfg, lim)
        save(fig, cfg['name'])
        plt.close(fig)
        print('%-18s %6.1f x %6.1f pt   scale ref %.2f pt (target %.2f)   '
              'box %s' % (cfg['name'], W, H, got, cfg['scale_ref'][2],
                          np.round(box, 3)))
    if COMBINE:
        try:
            from pypdf import PdfWriter, PdfReader, Transformation, PageObject
            pgs = [PdfReader(os.path.join(OUT, c['name'] + '.pdf')).pages[0]
                   for c in PANELS]
            gap = 14.0
            w = sum(float(p.mediabox.width) for p in pgs) + gap * (len(pgs) - 1)
            h = max(float(p.mediabox.height) for p in pgs)
            out = PageObject.create_blank_page(width=w, height=h)
            x = 0.0
            for p in pgs:
                out.merge_transformed_page(
                    p, Transformation().translate(x, h - float(p.mediabox.height)))
                x += float(p.mediabox.width) + gap
            wr = PdfWriter(); wr.add_page(out)
            with open(os.path.join(OUT, 'Figure2A_panelA.pdf'), 'wb') as fh:
                wr.write(fh)
            print('%-18s %6.1f x %6.1f pt' % ('Figure2A_panelA', w, h))
        except Exception as exc:
            print('combined page skipped:', exc)
    print('written to', OUT)


if __name__ == '__main__':
    main()
