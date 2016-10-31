LRA
===

Codes for Landcover Reconstruction Algorithm

The source code contains following functions:

REVEALS(pollen_counts, avg, alvc, u, Zmax, dwm)

Arguments

pollen_counts   contains EITHER data frame of single pollen count spreadsheet with following structure:
                - first row includes ages or depths or sample names
                - second row includes radius in metres
                - third row includes model type - 1 = Prentice bog model, 2 = Sugita lake model
                - fourth row and onwards include taxa with pollen counts
                OR list of data frames of multiple sites with the structure described above

avg             contains data frame where rows represent all calculated taxa, the first column is taxon name, the second column is alpha (pollen productivity estimate) and the third row is vg (fall speed of pollen)

alvc            contains data frame with variances and covariances of PPEs for all taxa arranged in columns and rows of the same size

u               wind speed in m/s
Zmax            maximum extent of the region in metres
dwm             type of dispersal model used: "gpm neutral", "lsm unstable"