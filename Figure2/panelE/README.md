# Figure 2E - receptor against pocket factor binding change

    python figure2_panelE.py

**Reads** `VDPV2_VP1_binding_energy_raw_data.xlsx`, sheet `Data`.

**Writes** into `figures/`: the panel as PDF, SVG and PNG, and
`marked_sets_positions.csv` with the coordinates of the highlighted sets.

Every Barcelona structural haplotype is one dot, its ddG against the CD155 receptor on
the horizontal axis and against the pocket factor on the vertical axis. Dashed
lines mark Sabin 2 at zero and dotted lines the +/- 1 REU criterion. The
fourteen sets that cross the criterion for both interactors are drawn on top of
the cloud, coloured where the pocket factor change also exceeds 2 REU and black
otherwise, each as its own group in the vector output so that labels and leader
lines can be added by hand.

**Needs** numpy, pandas, matplotlib, openpyxl.
