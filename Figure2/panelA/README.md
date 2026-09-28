# Figure 2A - haplotype-phenotype map, receptor against both antibodies

    python figure2_panelA.py

**Reads** `VDPV2_VP1_binding_energy_raw_data.xlsx`, sheet `Data`. The ddG of each
Barcelona structural haplotype against the CD155 receptor, mAb 9H2 and mAb 10D2, and the
haplotypes each set represents, are all taken from that one file.

**Writes** into `figures/`: an overview of all 808 structural haplotypes and a zoom on
the highlighted haplotypes, as PDF, SVG and PNG, plus both panels on one page.

Axes are 9H2 ddG, receptor ddG (negated, so the axis reads 3 to -3 as published)
and 10D2 ddG. Markers are small and semi-transparent so that overlapping points
in the dense part of the cloud stay distinguishable. Highlighted haplotypes are
drawn as a separate scatter on top of the cloud, so they arrive as their own
group in the vector output. Camera orientation and marker scale are fixed in the
CONFIG block.

**Needs** numpy, pandas, matplotlib, openpyxl, and pypdf to combine both panels.
