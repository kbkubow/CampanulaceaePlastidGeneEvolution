# CampanulaceaePlastidGeneEvolution
Data and R analysis script for manuscript "Identifying candidate genes for plastid-nuclear incompatibility: patterns of plastid gene evolution across the Campanulaceae"

Alignments folder contains the individual gene/gene family alignments used for analysis in PAML

concat_rooted.newick is the rooted version of the concatenated tree used as the input tree for PAML analyses and used for optimized branch length analyses using RaxML-NG

conat_unroot.newick is the unrooted concatenated tree used as the input tree for PAML analyses

dNdSdatasheet.csv is a comma separated data file containing the dN, dS, and dN/dS values extracted from PAML output, after dN/dS values were corrected for dN and dS or zero.

Analysis.R is the R script used to extract data from PAML output as well as visualize and analyze the dN/dS values. Contains code needed to reproduce the figures in the manuscript.
