LRA
===

## R codes for the Landscape Reconstruction Algorithm

This program calculates the REVEALS model [Sugita, 2007][1]. The calculation is possible for a single site or a combination of multiple sites. The function enables using deposition models for bogs and lakes, and two dispersal models (Gaussian plume model, Lagrangian stochastic model taken from the DISQOVER package [Theuerkauf,et al. 2016][2]).

~~~
Date (version)  : v.1.1, 4 Dec 2020
Author          : Petr Kuneš
Email           : petr.kunes@natur.cuni.cz
Citation        : Abraham, V., Oušková, V., & Kuneš, P. 2014. Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. PLoS ONE 9: e100117.
~~~
**Usage**

REVEALS(pollen_counts, avg, alvc, u, Zmax, dwm)

**Arguments**

<u>pollen_counts</u> – contains EITHER data frame of a single pollen count spreadsheet with the following structure:

* first row includes ages or depths or sample names (while using multiple sites colnames have to match at all sites)
* second row includes radius in metres
* third row includes model type: 1 = Prentice bog model, 2 = Sugita lake model
* fourth row and onwards include taxa with pollen counts

OR list of data frames of multiple sites with the structure described above

<u>avg</u> – contains data frame where rows represent all calculated taxa, the first column is taxon name, the second column is alpha (pollen productivity estimate) and the third row is vg (fall speed of pollen)

<u>alvc</u> – contains data frame with variances and covariances of PPEs for all taxa arranged in columns and rows of the same size

<u>u</u> – wind speed in m/s
<u>Zmax</u> – maximum extent of the region in metres
<u>dwm</u> – type of dispersal model used: "gpm neutral", "lsm unstable"

**Value**

Returns a list containing the following objects:

V_sites – vegetation estimates for particular sites

SE_sites – standard errors of vegetation estimantes for particular sites

Mean_V – mean vegetation estimates for all sites

Mean_SE – standard errors of mean vegetation estimates for all sites

No_sites – number of sites used for estimating particular time window



Please refer to 'example.R' for getting instructions how to run the function with the set of training data provided in 'example' folder.

## References:

[1]: Sugita, S. 2007. Theory of quantitative reconstruction of vegetation I: pollen from large sites REVEALS regional vegetation composition. *Holocene* 17: 229–241.

[2]: Theuerkauf, M., Couwenberg, J., Kuparinen, A., & Liebscher, V. 2016. A matter of dispersal: REVEALSinR introduces state-of-the-art dispersal models to quantitative vegetation reconstruction. *Vegetation History and Archaeobotany* 25: 541–553.