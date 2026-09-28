"""
Figure S1B - phylogeny of the most abundant VP1 amino acid haplotypes.

Aligns the sequences with MAFFT and infers the tree with IQ-TREE, letting
ModelFinder choose the substitution model and running 1000 ultrafast bootstrap
replicates.  The resulting Newick file is what the figure was drawn
from; the drawing itself was done in FigTree.

Run:      python figureS1_panelB.py
Data:     VP1_top58_haplotypes.fasta
Output:   VP1_top58_haplotypes_aligned.fasta and VP1_top58_haplotypes.treefile,
          beside the IQ-TREE report files
Requires: python >= 3.9, and mafft and iqtree2 (or iqtree) on the PATH

An alignment and a treefile from a previous run are already in this folder, so
the tree can be opened in FigTree without re-running anything.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

INPUT = 'VP1_top58_haplotypes.fasta'
ALIGNED = 'VP1_top58_haplotypes_aligned.fasta'
PREFIX = 'VP1_top58_haplotypes'

MAFFT_ARGS = ['--localpair', '--maxiterate', '1000']
IQTREE_ARGS = ['-m', 'MFP',        # ModelFinder picks the substitution model
               '-T', 'AUTO',       # threads
               '-B', '1000',       # ultrafast bootstrap replicates
               '--redo']


def need(binaries):
    for names in binaries:
        found = next((shutil.which(n) for n in names if shutil.which(n)), None)
        if found is None:
            sys.exit('%s not found on the PATH' % ' or '.join(names))
        yield found


def main():
    os.chdir(HERE)
    if not os.path.exists(INPUT):
        sys.exit('%s not found next to the script' % INPUT)
    mafft, iqtree = need([('mafft',), ('iqtree2', 'iqtree')])

    print('aligning with MAFFT ...')
    with open(ALIGNED, 'w') as out:
        subprocess.run([mafft] + MAFFT_ARGS + [INPUT], stdout=out, check=True)

    print('inferring the tree with IQ-TREE ...')
    subprocess.run([iqtree, '-s', ALIGNED, '-pre', PREFIX] + IQTREE_ARGS, check=True)

    tree = PREFIX + '.treefile'
    if not os.path.exists(tree):
        tree = PREFIX + '.contree'
    print('\ntree written to %s' % tree)
    print('open it in FigTree to reproduce the published panel')


if __name__ == '__main__':
    main()
