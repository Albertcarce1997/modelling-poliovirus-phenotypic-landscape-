"""
Figure 1A - frequency distribution of VP1 amino acid haplotypes.

Each haplotype is weighted by the number of reads supporting it, given by
size= in its header.  Haplotypes below 0.1 % of the reads are pooled into a
single grey slice, drawn last.  Slices are laid out clockwise from twelve
o'clock, largest first.

Run:      python figure1_panelA.py
Data:     VP1_amino_acid_haplotypes.fasta
Output:   figure1_panelA.png and figure1_panelA.csv
Requires: python >= 3.9, plotly, kaleido   (kaleido only to write the image)
"""
import csv
import os
import re
import sys

FASTA = 'VP1_amino_acid_haplotypes.fasta'
EXCLUDE = ['HT-612']              # excluded from the analysis, see the README
MIN_PERCENT = 0.1                 # below this a haplotype joins the pooled slice
POOLED_LABEL = 'Low frequency (<0.1%)'
POOLED_COLOUR = 'rgb(219, 217, 211)'
OUT_IMAGE = 'figure1_panelA.png'
OUT_TABLE = 'figure1_panelA.csv'
IMAGE_SIZE = 2400                 # pixels, square


HERE = os.path.dirname(os.path.abspath(__file__))


def read_counts(path):
    counts = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            size = re.search(r'size=(\d+)', line)
            name = re.search(r'(HT-\d+)', line)
            if size and name:
                counts[name.group(1)] = counts.get(name.group(1), 0) + int(size.group(1))
    return counts


def main():
    os.chdir(HERE)
    if not os.path.exists(FASTA):
        sys.exit('%s not found; run this script from the folder that contains it' % FASTA)
    counts = read_counts(FASTA)
    n_all = len(counts)
    for h in EXCLUDE:
        counts.pop(h, None)

    total = sum(counts.values())
    rows = sorted(({'Haplotype': h, 'Reads': n, 'Percent': 100.0 * n / total}
                   for h, n in counts.items()), key=lambda r: -r['Reads'])

    with open(OUT_TABLE, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['Haplotype', 'Reads', 'Percent'])
        w.writeheader()
        for r in rows:
            w.writerow({'Haplotype': r['Haplotype'], 'Reads': r['Reads'],
                        'Percent': round(r['Percent'], 4)})

    kept = [r for r in rows if r['Percent'] >= MIN_PERCENT]
    pooled = [r for r in rows if r['Percent'] < MIN_PERCENT]
    labels = [r['Haplotype'] for r in kept]
    values = [r['Reads'] for r in kept]
    if pooled:
        labels.append(POOLED_LABEL)
        values.append(sum(r['Reads'] for r in pooled))

    print('haplotypes in the file   %d' % n_all)
    print('excluded                 %s' % (', '.join(EXCLUDE) if EXCLUDE else 'none'))
    print('haplotypes plotted       %d' % len(rows))
    print('reads                    %d' % total)
    print('slices                   %d named + 1 pooled (%d haplotypes)'
          % (len(kept), len(pooled)))
    for r in kept[:5]:
        print('   %-10s %8d  %6.2f %%' % (r['Haplotype'], r['Reads'], r['Percent']))
    if pooled:
        print('   %-10s %8d  %6.2f %%' % ('pooled', sum(r['Reads'] for r in pooled),
                                          sum(r['Percent'] for r in pooled)))

    try:
        import plotly.graph_objects as go
        import plotly.colors as pc
    except ImportError:
        sys.exit('\nplotly is not installed; %s holds the values the pie shows' % OUT_TABLE)

    palette = pc.qualitative.Pastel + pc.qualitative.Pastel1 + pc.qualitative.Pastel2
    if len(kept) > len(palette):
        palette = palette + pc.qualitative.Light24
    colours = [palette[i % len(palette)] for i in range(len(kept))]
    if pooled:
        colours.append(POOLED_COLOUR)

    fig = go.Figure(data=[go.Pie(
        labels=labels, values=values, sort=False, direction='clockwise', rotation=0,
        textinfo='none', marker=dict(colors=colours))])
    fig.update_layout(showlegend=False, margin=dict(t=0, b=0, l=0, r=0),
                      paper_bgcolor='white')
    try:
        fig.write_image(OUT_IMAGE, width=IMAGE_SIZE, height=IMAGE_SIZE)
        print('\nwritten to %s and %s' % (OUT_IMAGE, OUT_TABLE))
    except Exception as exc:
        fig.write_html(OUT_IMAGE.replace('.png', '.html'))
        print('\nkaleido is not available (%s); wrote an HTML version instead' % exc)


if __name__ == '__main__':
    main()
