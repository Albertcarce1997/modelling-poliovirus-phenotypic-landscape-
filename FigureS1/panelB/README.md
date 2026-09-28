# Figure S1B - phylogeny of the most abundant haplotypes

    python figureS1_panelB.py

**Reads** `VP1_top58_haplotypes.fasta`: the 58 VP1 amino acid haplotypes at 0.1 %
of the reads or more, together with the NIE-ZAS-1 and Sabin 2 references.

**Writes** `VP1_top58_haplotypes_aligned.fasta` and
`VP1_top58_haplotypes.treefile`, beside the IQ-TREE report files.

The sequences are aligned with MAFFT (`--localpair --maxiterate 1000`) and the
tree inferred with IQ-TREE, with ModelFinder choosing the substitution model
(`-m MFP`) and 1000 ultrafast bootstrap replicates (`-B 1000`). The published
image was drawn from the resulting treefile in FigTree v1.4.5.

An alignment and a treefile from a previous run are already here, so the tree can
be opened in FigTree without re-running anything.
`VP1_top58_haplotype_attributes.txt` lists the read count of each haplotype, used
to annotate the tips.

**Needs** mafft and iqtree2 on the PATH, and FigTree to draw the tree.
