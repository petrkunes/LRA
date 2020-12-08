LRA
===

## R package for the Landscape Reconstruction Algorithm

This program calculates the REVEALS model [Sugita, 2007][2]. The calculation is possible for a single site or a combination of multiple sites. The function enables using deposition models for peatlands and lakes, and two dispersal models (Gaussian plume model [Prentice, 1985][1], Lagrangian stochastic model taken from the DISQOVER package [Theuerkauf,et al. 2016][3]).

~~~
Date (version)  : v0.1.0, 8 Dec 2020
Author          : Petr Kuneš
Email           : petr.kunes@natur.cuni.cz
Citation        : Abraham, V., Oušková, V., & Kuneš, P. 2014. Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. PLoS ONE 9: e100117.
~~~

How do you install the package from GitHub?

1. First, you need to install the `devtools` package. Run R and then type

   ```
   install.packages("devtools")
   ```

2. Load the `devtools` package.

   ```
   library(devtools)
   ```

3. In most cases, you just use `install_github()` function.

   ```
   install_github("petrkunes/LRA")
   ```

4. Load the `LRA` package

   ```
   library(LRA)
   ```



## References:

[1]: Prentice, I.C. 1985. Pollen representation, source area, and basin size: Toward a unified theory of pollen analysis. *Quaternary Research* 23: 76–86.

[2]: Sugita, S. 2007. Theory of quantitative reconstruction of vegetation I: pollen from large sites REVEALS regional vegetation composition. *Holocene* 17: 229–241.

[3]: Theuerkauf, M., Couwenberg, J., Kuparinen, A., & Liebscher, V. 2016. A matter of dispersal: REVEALSinR introduces state-of-the-art dispersal models to quantitative vegetation reconstruction. *Vegetation History and Archaeobotany* 25: 541–553.
