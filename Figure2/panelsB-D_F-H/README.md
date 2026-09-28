# Figures 2B to 2D and 2F to 2H - interface differences from Sabin 2

    python figure2_panelsB-D_F-H.py

**Reads** `Figure2_interface_descriptors.xlsx`: the buried surface area, the
number of contacting residues and the number of hydrogen bonds measured on every
retained model of each flex-ddG ensemble, one row per model. Sheet `Read_me`
documents the columns. Group `Sabin 2` is the reference the others are measured
from.

**Writes** into `figures/`: one grid per part of the figure, one row per measure,
each panel separately, and `plotted_differences.csv` with every plotted value.

Each haplotype is a dot at its mean difference from Sabin 2, with a bar from
mean - SD to mean + SD across the models of its ensemble; orange is above the
reference and blue below. Sabin 2 is the top row of each panel, in grey and at
zero by definition. Where a measure takes the same value in every model the bar
is absent and only the dot remains. All panels of a measure share one symmetric
axis; set `SCALE = 'panel'` for one axis per panel.

**Needs** numpy, pandas, matplotlib, openpyxl.
