# GMD-JGR-FESDIA
This is a published version of FESDIA. The archived version of FESDIA is hosted in zenodo as DOI. See [zenodo page](https://doi.org/10.5281/zenodo.15479943) for reference. 

## Installation
### System requirements

- R (version >= 4.0.0)
- Rtools (version >= 4.0.0)

Before installing the FESDIA, you will need to first install the following R-packages and their dependencies:

* *deSolve*, for the dynamic models developed in dtR;
* *FME*, to fit models to data, for sensitivity analysis, monte carlo simulations;
* *plot3D*, required for simple plotting


The code can be installed using: 
```r
devtools::install_github("stanleesocca/FESDIA")
```

And you can test the code with

```r
library(FESDIA)
sol <- FESDIAsolve()
plot(sol, which = c("FDET", "O2", "NO3", "SO4"), grid = FESDIAdepth())
```

Please cite the accompanying paper if the code is used in your work:

```
Nmor, Stanley Ifeanyi, Eric Viollier, Lucie Pastor, Bruno Lansard, Christophe Rabouille, and Karline Soetaert. "FESDIA (v1. 0): exploring temporal variations of sediment biogeochemistry under the influence of flood events using numerical modelling." Geoscientific Model Development Discussions 2022 (2022): 1-42.

Nmor, Stanley, Eric Viollier, Lucie Pastor, Bruno Lansard, and Christophe Rabouille. "Biogeochemical implication of massive episodic flood deposition: Model‐data integration." Journal of Geophysical Research: Oceans 130, no. 6 (2025): e2025JC022414.
```