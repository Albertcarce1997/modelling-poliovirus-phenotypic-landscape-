# Figure S1A - directional nucleotide substitution frequencies

    python figureS1_panelA.py

**Reads** `VP1_haplotypes_aligned_to_Sabin2.fasta`, the VP1 nucleotide haplotypes
aligned to the Sabin 2 reference. Each header carries `size=N`, the number of
reads supporting that haplotype.

**Writes** `substitution_frequencies.csv`: the twelve reference-to-haplotype
changes with their weighted counts and their percentage of all substitutions.

Every haplotype is compared with the reference column by column, over the span
the reference covers. Each difference counts once per supporting read, so a
substitution carried by an abundant haplotype weighs more than the same
substitution in a rare one. Transitions are A<->G and C<->T; everything else is
a transversion.

The alignment covers 1,210 haplotypes and gives a transition/transversion ratio
of 14.09. The published panel draws these percentages as arrows between the four
bases.

**Needs** nothing beyond the Python standard library.
