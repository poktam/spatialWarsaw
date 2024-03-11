# spatialWarsaw R package
Repository for spatialWarsaw project in R

The functions collected in the spatialWarsaw package support spatial analysis on geolocalised points. They use spatial machine learning, spatial statistics and spatial econometric approaches. You can detect density clusters (with `QDC`, `ClustCont`, `ClustDisjoint`, `bootdbscan` functions), measure global and local agglomeration of points (`SPAG`, `ETA`, `FLE` functions), construct the spatial weight matrix W from Voronoi polygons (`tessW` function), determine the best knn structure of W with AIC (`bestW` function), check the correlation between spatial lags with different knn (`corrSpatialLags` function), and study semi-variance by expanding knn (`semiVarKnn` function). You can test and generate spatial point patterns that follow Benford's law (with `SpatBenfordTest`, `SpatBenfordPattern`). In spatial econometrics, it can run switching regime models - regression in density subgroups (`ssr` function) and bootstrapped spatial regression (`BootSpatReg` function), generate out-of-sample predictions (`SpatPredTess` function), approximate standard errors for large samples (`ApproxSERoot2` function), and rasterise and cluster GWR coefficients (with `rastClustGWR` function) to check their spatiotemporal stability (`STS` function).  
***Please note that the package is still under development, although all the functions provided work. Currently we are working on cleaning up the code and creating help for specific functions. The package will soon be updated with sample data.***

## Installation

To install the current **development** version directly from Github (this repository), the following commands can be used:

``` r
# install.packages("devtools")
devtools::install_github ("poktam/spatialWarsaw")
```
or

``` r
# install.packages("remotes")
remotes::install_github ("poktam/spatialWarsaw")
```

To load the package and check the version:

``` r
library (spatialWarsaw)
packageVersion ("spatialWarsaw")
```

## Authors and Contributors
The package was created thanks to the work of former and current members and collaborators of the [Spatial Warsaw Team](https://spatial.wne.uw.edu.pl/) affiliated to the Faculty of Economic Sciences of the University of Warsaw.

- Katarzyna Kopczewska
- Mateusz Kopyt
- Maria Kubara
- Ewa Dobrowolska

## License

The spatialWarsaw package as free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as published by 
the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; 
without even the implied warranty of merchantability or fitness for a particular purpose.
See the [GPL-3 (GNU General Public License version 3)](https://www.gnu.org/licenses/gpl-3.0.en.html)
for more details.

A copy of the GNU General Public License, version 3, is available also in this package repository and at <https://www.r-project.org/Licenses/GPL-3>.
