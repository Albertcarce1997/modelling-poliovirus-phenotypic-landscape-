# Figures 1C, 1D and 1E - VP1 structure

Rendered with UCSF ChimeraX.

    chimerax --script "figure1_panelsC-E.py --input_file AA_diversity.xlsx"

**Reads** `AA_diversity.xlsx`, per-position amino acid diversity for the 301 VP1
positions, and `VP1_Sabin_cristal_cleaned.pdb`, the Sabin 2 VP1 structure.

**Writes** a coloured PDB and the ChimeraX views. Panel C maps diversity onto the
structure, darker meaning more diverse. Panels D and E highlight the four
antigenic regions and the receptor-binding residues, using the residue lists in
Table S1 of the manuscript.

**Needs** ChimeraX. The script runs inside it; no Python environment is required.
