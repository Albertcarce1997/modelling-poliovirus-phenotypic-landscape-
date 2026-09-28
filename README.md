# VDPV2 VP1 figures - data and code

Data and code for the figures of *Mutational and phenotypic landscape of
circulating Vaccine-Derived Poliovirus 2 in sewage in a high vaccine coverage
area*.

Each folder is self-contained: it holds one script, the data that script reads
and a README describing both. Download a folder and run the script; it finds its
data beside itself and writes its output there, whatever the working directory.
No folder depends on any other. Scripts are named after the panels they draw.

| Folder | Panels | Tool |
|---|---|---|
| `Figure1/panelA` | 1A | Python, Plotly |
| `Figure1/panelB` | 1B | GraphPad Prism |
| `Figure1/panelsC-E` | 1C, 1D, 1E | UCSF ChimeraX |
| `Figure1/panelsF-Q` | 1F-1Q dot plots | Python, matplotlib |
| `Figure1/stacked_bars` | bars inside 1F-1Q | GraphPad Prism |
| `Figure2/panelA` | 2A | Python, matplotlib |
| `Figure2/panelsB-D_F-H` | 2B-2D, 2F-2H | Python, matplotlib |
| `Figure2/panelE` | 2E | Python, matplotlib |
| `FigureS1/panelA` | S1A | Python |
| `FigureS1/panelB` | S1B | MAFFT, IQ-TREE, FigTree |
| `FigureS2` | S2 | Python |
| `FigureS3` | S3 | Python |

Published panels were assembled and annotated in Adobe Illustrator, so colours,
labels and leader lines may differ from the raw script output. The values are
the same.

## Requirements

Python 3.9 or later. Most folders need only `numpy`, `pandas`, `matplotlib` and
`openpyxl`; `Figure1/panelA` needs `plotly` and `kaleido`, and
`FigureS1/panelA` needs nothing beyond the standard library. Each README states
what its own script needs.

    pip install numpy pandas matplotlib openpyxl plotly kaleido

The matplotlib figures use Arial; any sans-serif font works, and the fallback
list is the first setting in each script's CONFIG block.
