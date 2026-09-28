# Stacked-bar plots inside Figures 1F to 1Q

Plotted with GraphPad Prism. One file per dataset:

    Stacked_Bar_plots_sewage.prism    Barcelona sewage
    Stacked_Bar_plots_Nigeria.prism   Nigeria cVDPV2
    Stacked_Bar_plots_iVDPV.prism     iVDPV2

`stacked_bar_values.csv` holds the plotted values for all three datasets, so the
bars can be read without Prism. For each dataset, interactor, VP1 site and
direction of the effect it gives:

- `Structural_haplotypes_with_a_replacement_at_the_site`, the numerator
- `Structural_haplotypes_crossing_the_threshold`, the denominator for that
  interactor, that is, how many structural haplotypes reach |ddG| > 1 REU
- `Percent`, the height of the bar

Sites are the four antigenic sites,
the receptor site and the pocket factor site, as listed in Table S1 of the
manuscript.
