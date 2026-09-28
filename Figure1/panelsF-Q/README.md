# Figures 1F to 1Q - predicted changes in VP1 binding energy

    python figure1_panelsF-Q.py

**Reads** `VDPV2_VP1_binding_energy_raw_data.xlsx`. Sheet `Read_me` documents
every column; sheet `Data` holds one row per VP1 structural haplotype and interactor,
with the individual Rosetta flex-ddG trajectory energies and everything derived
from them.

**Writes** into `figures/`: one four-panel figure per dataset, each panel
separately, and `marked_sets_positions.csv` with the rank and ddG of the
labelled sets. PDF, SVG and PNG, with text kept as text.

Each structural haplotype is one dot, sorted by ddG up the vertical axis, so that axis is
rank and the horizontal axis is ddG. All twelve panels share one symmetric axis
from -15 to +15 REU. Dashed lines mark the +/- 1 REU criterion; dots are dark
green below -1, grey within and red beyond +1, and the counts in the upper left
give the number in each group.

The stacked-bar plots embedded in these panels are in `Figure1/stacked_bars`.

**Needs** numpy, pandas, matplotlib, openpyxl.
