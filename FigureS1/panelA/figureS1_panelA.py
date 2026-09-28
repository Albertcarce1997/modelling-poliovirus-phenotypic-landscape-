"""
Figure S1, panel A - directional nucleotide substitution frequencies.

Every VP1 nucleotide haplotype is compared with the Sabin 2 reference, column by column,
over the span of the alignment covered by the reference.  Each difference is
counted once per read supporting the haplotype, so a substitution carried by an
abundant haplotype weighs more than the same substitution in a rare one.  The
twelve reference-to-haplotype changes are reported as a percentage of all
substitutions, split into transitions and transversions.

Run:      python figureS1_panelA.py
Data:     VP1_haplotypes_aligned_to_Sabin2.fasta
Output:   substitution_frequencies.csv
Requires: python >= 3.9 (no third-party packages)
"""
import csv
import os
import re
import sys

FASTA = 'VP1_haplotypes_aligned_to_Sabin2.fasta'
REFERENCE = 'Sabin_2_VP1'
OUT = 'substitution_frequencies.csv'

BASES = ('A', 'C', 'G', 'T')
PURINES = {'A', 'G'}


HERE = os.path.dirname(os.path.abspath(__file__))


def read_fasta(path):
    records, header, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks)))
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, ''.join(chunks)))
    return records


def read_count(header):
    """Number of reads supporting a haplotype, from size=N in its header."""
    m = re.search(r'size=(\d+)', header)
    return int(m.group(1)) if m else 1


def kind(a, b):
    return 'transition' if (a in PURINES) == (b in PURINES) else 'transversion'


def main():
    os.chdir(HERE)
    if not os.path.exists(FASTA):
        sys.exit('%s not found; run this script from the folder that contains it' % FASTA)
    records = read_fasta(FASTA)
    ref = next((s for h, s in records if h.split()[0] == REFERENCE), None)
    if ref is None:
        sys.exit('reference %s not found in %s' % (REFERENCE, FASTA))
    ref = ref.upper()

    # restrict to the span the reference actually covers
    first = min(i for i, c in enumerate(ref) if c != '-')
    last = max(i for i, c in enumerate(ref) if c != '-')

    counts = {(a, b): 0 for a in BASES for b in BASES if a != b}
    reads = 0
    for header, seq in records:
        if header.split()[0] == REFERENCE:
            continue
        seq = seq.upper()
        n = read_count(header)
        reads += n
        for i in range(first, min(last + 1, len(seq))):
            a, b = ref[i], seq[i]
            if a not in BASES or b not in BASES or a == b:
                continue
            counts[(a, b)] += n

    total = sum(counts.values())
    ts = sum(v for (a, b), v in counts.items() if kind(a, b) == 'transition')
    tv = total - ts

    rows = []
    for a in BASES:
        for b in BASES:
            if a == b:
                continue
            v = counts[(a, b)]
            rows.append({'From': a, 'To': b, 'Type': kind(a, b), 'Substitutions': v,
                         'Percent_of_all_substitutions': round(100.0 * v / total, 4) if total else 0.0})
    rows.sort(key=lambda r: -r['Substitutions'])

    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['From', 'To', 'Type', 'Substitutions',
                                           'Percent_of_all_substitutions'])
        w.writeheader()
        w.writerows(rows)

    print('haplotypes compared      %d' % (len(records) - 1))
    print('reads represented        %d' % reads)
    print('substitutions (weighted) %d' % total)
    print('transitions              %d (%.2f%%)' % (ts, 100.0 * ts / total))
    print('transversions            %d (%.2f%%)' % (tv, 100.0 * tv / total))
    print('transition/transversion  %.2f' % (ts / tv if tv else float('inf')))
    print()
    for r in rows:
        print('  %s>%s  %-12s %12d  %6.2f %%' % (r['From'], r['To'], r['Type'],
                                                 r['Substitutions'],
                                                 r['Percent_of_all_substitutions']))
    print('\nwritten to %s' % OUT)


if __name__ == '__main__':
    main()
