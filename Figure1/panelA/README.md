# Figure 1A - haplotype frequency distribution

    python figure1_panelA.py

**Reads** `VP1_amino_acid_haplotypes.fasta`, the 1,030 VP1 amino acid haplotypes
of the Barcelona sewage quasispecies. Each header carries `size=N`, the number of
reads supporting that amino acid haplotype.

**Writes** `figure1_panelA.csv`, every amino acid haplotype with its read count and
percentage, and `figure1_panelA.png`, the pie chart.

HT-612 is present in the fasta but excluded from the figure, as it is from the
rest of the analysis: it carries several consecutive replacements consistent
with an artefactual alignment. `EXCLUDE` at the top of the script controls this.
The remaining 1,029 amino acid haplotypes account for 390,362 reads, of which 58
reach 0.1 % or more; the rest are pooled into one grey slice drawn last. Slices
run clockwise from twelve o'clock, largest first. Labels are not drawn, since the
published panel is annotated separately.

**Needs** `plotly` and `kaleido` for the image. Without them the script still
writes the CSV, which holds every value the pie shows.
