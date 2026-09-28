# Figure S2 - mean quality score along the VP1 gene

    python figureS2.py

**Reads** `position_qscores.csv`: one row per VP1 nucleotide position, with the
number of base calls covering it and the mean and standard deviation of their
Phred quality scores, pooled over every retained amino acid haplotype.

**Writes** `figureS2.pdf`, `.svg` and `.png`.

The 909 positions have an overall mean of Q40.7; the lowest, Q30.1, is at
position 244.

**Needs** pandas and matplotlib.
