LRA
===

Codes for Landcover Reconstruction Algorithm

The source code contains following functions:

REVEALS(file_name_list, file_name_avg, file_alvc, u, Zmax)

Arguments

file_name_list  contains EITHER name of comma delimited file with following structure:
                - first row includes ages or depths or sample names
                - second row includes radius in metres
                - third row includes model type - 1 = Prentice, 2 = Sugita
                - fourth row and onwards include taxa with pollen counts
                OR name of text file with a list of filenames of input files with the structure described above

file_name_avg   contains name of comma delimited file where rows represent all taxa calculated, the first column is taxon name,                 the second column is alpha (pollen productivity estimate) and the third row is vg (fall speed of pollen)

file_alvc       contains name of comma delimited file with variances and covariances of PPE's for all taxa arranged in columns                  and rows of the same size

u               wind speed in m/s
Zmax            maximum extent of the region in km
