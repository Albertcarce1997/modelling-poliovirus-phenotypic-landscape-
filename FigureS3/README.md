# Figure S3 - distance of the haplotypes from the most abundant one

    python figureS3.py

**Reads** `haplotype_mutation_counts.csv`: one row per nucleotide haplotype, with
its read count and the number of nucleotide differences separating it from HT-1,
the most abundant haplotype.

**Writes** `figureS3.pdf`, `.svg` and `.png`, and `figureS3_values.csv` with the
plotted values.

Bars give the number of haplotypes at each distance, with that count above them;
the line gives the share of all reads those haplotypes represent, on the right
axis. HT-1 itself is not drawn but stays in the read denominator; set
`DROP_ZERO = False` to include it.

The fasta files the counts were derived from, `VP1_nucleotide_haplotypes.fasta`
and `HT-1_reference.fasta`, are included for reference.

**Needs** pandas and matplotlib.
