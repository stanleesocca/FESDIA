# Functions to retrieve variables for the PHDIA and FESDIA model.

`FESDIA0D, FESDIA1D, FESDIAsvar` retrieve the (0-dimensional or
1-dimensional) output variables or state variables of FESDIA model
solutions. When called with `out` not specified, will return the names
of these variables and their units.

`PHDIA0D, PHDIA1D, PHDIAsvar` retrieve the (0-dimensional or
1-dimensional) output variables or state variables of PHDIA model
solutions. When called with `out` not specified, will return the names
of these variables and their units.

`FESDIAextract` extracts model output variables of any dimension.

`FESDIApH` estimates the pH profile(s) from FESDIA model solutions.

## Usage

``` r
FESDIA0D(out, as.vector = FALSE, which = NULL) 
  FESDIA1D(out, which = NULL) 
  FESDIAsvar(out, which = NULL) 
  
  PHDIA0D(out, as.vector = FALSE, which = NULL) 
  PHDIA1D(out, which = NULL) 
  PHDIAsvar(out, which = NULL) 
  
  FESDIAextract(out, which = NULL) 
  FESDIApH(out)
```

## Arguments

- out :

  an output object returned by
  [FESDIAsolve](https://stanleesocca.github.io/FESDIA/reference/FESDIAsolve.md)
  or
  [FESDIAdyna](https://stanleesocca.github.io/FESDIA/reference/FESDIAdyna.md).
  If `NULL`, `FESDOIAparms` will return the default (parameter) values.

- as.vector :

  if `TRUE` will return the parameter vector, else a data.frame that
  also contains the units.

- which :

  if not `NULL`, a vector with names of the variables to return.

## Value

`FESDIA0D` and `FESDIA1D` return the output variables of the solution as
a vector or data.frame. For dynamic runs, the output is averaged over
the run.

`FESDIA1D` always returns the sediment depth and the porosity as the
first two columns.

## Author

Karline Soetaert

## References

Soetaert K, PMJ Herman and JJ Middelburg, 1996a. A model of early
diagenetic processes from the shelf to abyssal depths. Geochimica
Cosmochimica Acta, 60(6):1019-1040.

Soetaert K, PMJ Herman and JJ Middelburg, 1996b. Dynamic response of
deep-sea sediments to seasonal variation: a model. Limnol. Oceanogr.
41(8): 1651-1668.

## Examples

``` r
# =====================
# defaults
# =====================

  FESDIAsvar()
#>     names             units                         description
#> 1    FDET    mmolC/m3 solid      Fast decaying Detritus (solid)
#> 2    SDET    mmolC/m3 solid      Slow decaying Detritus (solid)
#> 3      O2   mmolO/m3 liquid                     Oxygen (liquid)
#> 4     NO3   mmolN/m3 liquid                    Nitrate (liquid)
#> 5     NO2   mmolN/m3 liquid                    Nitrite (liquid)
#> 6     NH3   mmolN/m3 liquid           Ammonium/ammonia (liquid)
#> 7     DIC   mmolC/m3 liquid Dissolved Inorganic Carbon (liquid)
#> 8      Fe  mmolFe/m3 liquid                       Fe2+ (liquid)
#> 9   FeOH3   mmolFe/m3 solid                    Fe-oxide (solid)
#> 10    H2S   mmolS/m3 liquid                   Sulphide (liquid)
#> 11    SO4   mmolS/m3 liquid                   Sulphate (liquid)
#> 12    CH4 mmolCH4/m3 liquid                    Methane (liquid)
#> 13    PO4   mmolP/m3 liquid                  Phosphate (liquid)
#> 14    FeP    mmolP/m3 solid                Iron-bound P (solid)
#> 15    CaP    mmolP/m3 solid                  Ca-bound P (solid)
#> 16   Pads    mmolP/m3 solid                  Adsorbed P (solid)
#> 17    ALK mmolALK/m3 liquid                 Alkalinity (liquid)
#> 18 FeOH3B   mmolFe/m3 solid         Crystalline Fe-oxid (solid)
#> 19     Mn  mmolMn/m3 liquid                       Mn2+ (liquid)
#> 20   MnO2   mmolMn/m3 solid                    Mn-oxide (solid)
#> 21  MnO2B   mmolMn/m3 solid        Crystalline Mn-oxide (solid)
  head(FESDIA0D())
#>         names values        units               description
#> 1      O2flux     NA nmolO2/cm2/d  O2 influx sediment-water
#> 2  O2deepflux     NA nmolO2/cm2/d  O2 efflux lower boundary
#> 3     NO3flux     NA  nmolN/cm2/d NO3 influx sediment-water
#> 4 NO3deepflux     NA  nmolN/cm2/d NO3 efflux lower boundary
#> 5     NO2flux     NA  nmolN/cm2/d NO2 influx sediment-water
#> 6 NO2deepflux     NA  nmolN/cm2/d NO2 efflux lower boundary
  FESDIA1D()
#>              names               units
#> 1              TOC                   %
#> 2       DICprodMin  nmolC/cm3 liquid/d
#> 3       DINprodMin  nmolN/cm3 liquid/d
#> 4       DIPprodMin  nmolP/cm3 liquid/d
#> 5           O2prod  nmolO/cm3 liquid/d
#> 6          Oxicmin  nmolC/cm3 liquid/d
#> 7       Denitrific  nmolC/cm3 liquid/d
#> 8         Feredmin  nmolC/cm3 liquid/d
#> 9           BSRmin  nmolC/cm3 liquid/d
#> 10         Methmin  nmolC/cm3 liquid/d
#> 11          nitri1  nmolN/cm3 liquid/d
#> 12          nitri2  nmolN/cm3 liquid/d
#> 13         Anammox  nmolN/cm3 liquid/d
#> 14          Feoxid nmolFe/cm3 liquid/d
#> 15         H2Soxid  nmolS/cm3 liquid/d
#> 16         CH4oxid  nmolC/cm3 liquid/d
#> 17             AOM  nmolS/cm3 liquid/d
#> 18         FeSprod nmolFe/cm3 liquid/d
#> 19       FePadsorp nmolFe/cm3 liquid/d
#> 20       FePdesorp   nmolP/cm3 solid/d
#> 21         CaPprod  nmolP/cm3 liquid/d
#> 22         CaPdiss   nmolP/cm3 solid/d
#> 23         Padsorb   nmolP/cm3 solid/d
#> 24       H2Soxsurf  nmolS/cm3 liquid/d
#> 25       CH4oxsurf  nmolC/cm3 liquid/d
#> 26   O2distConsump  nmolO/cm3 liquid/d
#> 27         ALKprod   nmol/cm3 liquid/d
#> 28      DICprodCH4  nmolC/cm3 liquid/d
#> 29        MPBCprod   nmolC/cm3 solid/d
#> 30    MPBuptakeNO3  nmolN/cm3 liquid/d
#> 31    MPBuptakeNH3  nmolN/cm3 liquid/d
#> 32    MPBuptakePO4  nmolP/cm3 liquid/d
#> 33    MPBuptakeDIC  nmolC/cm3 liquid/d
#> 34        MnredMin nmolMn/cm3 liquid/d
#> 35         sumMnO2    nmolMn/cm3 solid
#> 36        sumFeOH3  nmolFe/cm3 solid/d
#> 37     FeOxidMnASC nmolFe/cm3 liquid/d
#> 38 H2SOxidFeOH3ASC  nmolS/cm3 liquid/d
#> 39  H2SoxidMnO2ASC nmolS/cm3 liquid /d
#> 40         MnSprod nmolMn/cm3 liquid/d
#> 41          MnOxid nmolMn/cm3 liquid/d
#> 42         AgeFeOx  nmolFe/cm3 solid/d
#> 43         AgeMnOx  nmolMn/cm3 solid/d
#> 44       MnCO3prec nmolMn/cm3 liquid/d
#>                                     description
#> 1                Total Organic Carbon % profile
#> 2       DIC production profile (mineralisation)
#> 3       DIN production profile (mineralisation)
#> 4       DIP production profile (mineralisation)
#> 5     O2 production profile (microphytobenthos)
#> 6                   Oxic mineralisation profile
#> 7                       Denitrification profile
#> 8           Fe reduction mineralisation profile
#> 9     Sulphate reduction mineralisation profile
#> 10         Methanogensis mineralisation profile
#> 11 Nitrification step 1 profile (NH3 oxidation)
#> 12 Nitrification step 2 profile (NO2 oxidation)
#> 13                              Anammox profile
#> 14                       Fe2+ oxidation profile
#> 15                         H2S oxidaton profile
#> 16                        CH4 oxidation profile
#> 17       Anaerobic oxidation of methane profile
#> 18                       FeS production profile
#> 19                       FeP adsorption profile
#> 20                       FeP desorption profile
#> 21                       CaP production profile
#> 22                      CaP dissolution profile
#> 23                         P adsorption profile
#> 24        H2S oxidation with surface O2 profile
#> 25        CH4 oxidation with surface O2 profile
#> 26  O2 uptake oxidation with surface O2 profile
#> 27                Alkalinity production profile
#> 28           DIC production via Methane profile
#> 29                       MPB production profile
#> 30                       MPB NO3 uptake profile
#> 31                       MPB NH3 uptake profile
#> 32                       MPB PO4 uptake profile
#> 33                       MPB DIC uptake profile
#> 34          Mn reduction mineralisation profile
#> 35                           Ascorbic Mangenese
#> 36                                Ascorbic Iron
#> 37                           Fe oxidation by Mn
#> 38                         H2S oxidation by Fe3
#> 39                         H2S oxidation by Mn4
#> 40                          MnS profile profile
#> 41                                 Mn oxidation
#> 42                          FeOx Ageing profile
#> 43                          MnOx Ageing profile
#> 44                          MnCO3 precipitation
  
# =====================
# variables of runs  
# =====================

# defaults
  defsteady <- FESDIAsolve()
  defdyn    <- FESDIAdyna()

# altered steady-state run
  out  <- FESDIAsolve(parms = list(Cflux = 10))
  dyna <- FESDIAdyna(parms = list(Cflux = 10), CfluxForc = list(amp = 1))

# 0-D outputs
  cbind(steady.default = FESDIA0D(defsteady, as.vector = TRUE), 
        dyna.default   = FESDIA0D(defdyn, as.vector = TRUE), 
        out            = FESDIA0D(out,  as.vector = TRUE),
        dyna           = FESDIA0D(dyna, as.vector = TRUE))
#>                 steady.default   dyna.default            out           dyna
#> O2flux            5.069599e+02   5.069599e+02   1.727286e+01   1.740987e+01
#> O2deepflux       1.041187e-161  2.624450e-171   3.864035e-04   3.864034e-04
#> NO3flux          -2.299858e+01  -2.299858e+01  -2.067288e+00  -2.078977e+00
#> NO3deepflux       2.820115e-07   2.820115e-07   1.825395e-05   1.825395e-05
#> NO2flux          -1.062772e+01  -1.062772e+01  -2.086511e+00  -2.088189e+00
#> NO2deepflux       4.678466e-67   4.678466e-67   0.000000e+00   2.376743e-35
#> NH3flux          -2.597988e+01  -2.597988e+01   2.658362e+00   2.654881e+00
#> NH3deepflux       1.223195e-04   1.223195e-04   0.000000e+00   5.354006e-35
#> PO4flux          -4.307725e+00  -4.307725e+00  -9.429387e-02  -9.545796e-02
#> PO4deepflux       2.015280e-05   2.015280e-05   1.108627e-06   1.108627e-06
#> DICflux          -4.550175e+02  -4.550175e+02  -9.996355e+00  -1.011635e+01
#> DICdeepflux       4.502289e-03   4.502289e-03   2.912817e-03   2.912817e-03
#> Feflux           -5.522569e-02  -5.522569e-02  -6.537165e-04  -6.564688e-04
#> Fedeepflux        1.403179e-29   8.327928e-34   0.000000e+00  -4.652492e-46
#> H2Sflux          -6.904489e+00  -6.904489e+00  -3.297232e-07  -3.308020e-07
#> H2Sdeepflux       2.451697e-04   2.451697e-04   0.000000e+00   9.734673e-45
#> SO4flux           9.611513e+00   9.611513e+00   4.246700e-02   4.246701e-02
#> SO4deepflux       4.198996e-02   4.198996e-02   4.246581e-02   4.246581e-02
#> CH4flux          -4.391394e-08  -4.391394e-08  -9.365122e-10  -9.390685e-10
#> CH4deepflux       1.702274e-44   1.702274e-44   0.000000e+00  -5.390752e-48
#> ALKflux          -7.470144e+00  -7.470144e+00   6.910377e+00   6.922399e+00
#> ALKdeepflux       4.262150e-03   4.262150e-03   3.277133e-03   3.277133e-03
#> FDETflux          4.109589e+02   4.109589e+02   9.000000e+00   9.000000e+00
#> FDETdeepflux     1.001633e-191  1.001633e-191  2.193575e-193  2.193575e-193
#> SDETflux          4.566210e+01   4.566210e+01   1.000000e+00   1.000000e+00
#> SDETdeepflux      1.168446e-71   1.168446e-71   2.558897e-73   2.558897e-73
#> FePsurfflux       0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> FePdeepflux       1.367821e-10   3.597467e-10   4.463963e-05   4.463963e-05
#> CaPsurfflux       0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> CaPdeepflux       0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> FeOH3surfflux     5.000000e-01   5.000000e-01   5.000000e-01   5.000000e-01
#> FeOH3deepflux     2.949206e-55  1.574236e-180   0.000000e+00   1.289351e-36
#> OrgCflux          4.566210e+02   4.566210e+02   1.000000e+01   1.000000e+01
#> OrgNflux          6.892393e+01   6.892393e+01   1.509434e+00   1.509434e+00
#> OrgPflux          4.307745e+00   4.307745e+00   9.433962e-02   9.433962e-02
#> DINDIPflux        1.136991e+01   1.136991e+01  -6.268416e+00  -2.411419e+01
#> DINDIPmean        2.654843e+00   2.654843e+00   1.654154e+01   1.653268e+01
#> DINDIPdeep        2.652951e+00   2.652951e+00   1.646536e+01   1.646536e+01
#> TotMin            4.566210e+02   4.566210e+02   1.000000e+01   1.009290e+01
#> TotOxic           4.270496e+02   4.270496e+02   9.987815e+00   1.008053e+01
#> TotDenit          9.437490e+00   9.437490e+00   9.012522e-03   9.172229e-03
#> TotFered          1.253984e-01   1.253984e-01   8.228437e-04   8.304999e-04
#> TotBSR            1.927247e+01   1.927247e+01   2.395343e-06   2.409533e-06
#> TotMeth           6.393563e-01   6.393563e-01   7.871787e-08   7.918422e-08
#> PartOxic          9.352385e-01   9.352385e-01   9.987815e-01   9.987981e-01
#> PartDenit         2.066810e-02   2.066810e-02   9.012522e-04   9.166109e-04
#> PartFered         2.746225e-04   2.746225e-04   8.228437e-05   7.498517e-05
#> PartBSR           4.220671e-02   4.220671e-02   2.395343e-07   2.399077e-07
#> PartMethano       1.400190e-03   1.400190e-03   7.871787e-09   7.884058e-09
#> TotNitri1         4.206010e+01   4.206010e+01   4.164412e+00   4.174930e+00
#> TotNitri2         3.054857e+01   3.054857e+01   2.074517e+00   2.083290e+00
#> TotAnammox        8.838149e-01   8.838149e-01   3.383928e-03   3.474786e-03
#> TotFeoxid         3.031691e+00   3.031691e+00   2.556193e-03   2.582738e-03
#> TotH2Soxid        3.684392e-01   3.684392e-01   1.193643e-09   1.196895e-09
#> TotCH4oxid        1.795032e-02   1.795032e-02   3.841795e-08   3.864853e-08
#> TotAOM            3.017278e-01   3.017278e-01   4.474724e-12   4.514418e-12
#> TotFeSprod        6.250444e-02   6.250444e-02   5.662776e-16   8.425100e-16
#> TotFePprod        3.280683e-03   3.280683e-03   2.713154e-04   2.726046e-04
#> TotCaPprod        0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotFePdesorp      2.827091e-03   2.827084e-03   1.134290e-04   1.148430e-04
#> TotCaPdiss        0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotPadsorb        0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotNH3prod        6.892393e+01   6.892393e+01   1.509434e+00   1.523456e+00
#> TotPO4prod        4.307745e+00   4.307745e+00   9.433962e-02   9.521601e-02
#> TotNH3ads         1.468435e+01   1.468435e+01  -1.502552e+00  -1.500623e+00
#> TotO2prod         0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotH2Soxsurf      0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotCH4oxsurf      0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotALkprod        7.474406e+00   7.474406e+00  -6.907099e+00  -6.914882e+00
#> PartPremoved      1.052969e-04   1.052987e-04   1.673596e-03   4.240546e-03
#> PartNremoved      1.351870e-01   1.351870e-01   9.260341e-03   1.223942e-02
#> TotMPBNO3uptake   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotMPBNH3uptake   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotMPBPO4uptake   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotMPBDICuptake   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotMPBO2prod      0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotalFDET         6.000000e+03   6.000000e+03   1.314000e+02   1.326430e+02
#> TotalSDET         3.333333e+05   3.333333e+05   7.300000e+03   7.356659e+03
#> TotalO2           1.018914e+02   1.018914e+02   1.415583e+04   1.415051e+04
#> TotalNO3          4.914866e+01   4.914866e+01   6.650106e+02   6.656896e+02
#> TotalNO2          2.094459e+00   2.094459e+00   1.040778e-01   1.045185e-01
#> TotalNH3          1.864803e+03   1.864803e+03   2.089222e-01   2.094506e-01
#> TotalDIC          1.634170e+05   1.634170e+05   1.065439e+05   1.065498e+05
#> TotalFe           1.038276e+00   1.038276e+00   2.874307e-05   2.913416e-05
#> TotalFeOH3        7.688476e+02   7.688476e+02   3.638801e+07   3.638801e+07
#> TotalH2S          8.659309e+03   8.659309e+03   8.068827e-09   8.115937e-09
#> TotalSO4          1.536889e+06   1.536889e+06   1.553720e+06   1.553720e+06
#> TotalCH4          3.274970e-01   3.274970e-01   4.811529e-12   4.854210e-12
#> TotalPO4          7.217172e+02   7.217172e+02   4.022138e+01   4.028417e+01
#> TotalFeP          2.023936e+00   2.031286e+00   1.625477e+03   1.625467e+03
#> TotalCaP          0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> TotalPads         0.000000e+00  1.699419e-158   0.000000e+00   0.000000e+00
#> FeOH3Bsurfflux    5.000000e-01   5.000000e-01   5.000000e-01   5.000000e-01
#> FeOH3Bdeepflux    5.635807e-56  5.685062e-179   9.993281e-01   9.993281e-01
#> Mnflux           -2.829347e-01  -2.829347e-01  -3.727186e-03  -3.746189e-03
#> Mndeepflux        1.520882e-62   1.520882e-62   0.000000e+00   5.310044e-41
#> MnO2surfflux      5.000000e-01   5.000000e-01   5.000000e-01   5.000000e-01
#> MnO2deepflux     2.378218e-178  2.378218e-178   0.000000e+00   2.329690e-39
#> MnO2Bsurfflux     5.000000e-01   5.000000e-01   5.000000e-01   5.000000e-01
#> MnO2Bdeepflux    1.869506e-172  1.869506e-172   9.955591e-01   9.955591e-01
#> TotMnred          9.673114e-02   9.673114e-02   2.347379e-03   2.362180e-03
#> TotMnOxid         3.036029e-02   3.036029e-02   2.864152e-04   2.886147e-04
#> TotMnSprod        3.379908e-04   3.379908e-04   1.171191e-16   1.718486e-16
#> TotMnASC          3.039111e+02   3.039111e+02   3.625077e+07   3.625077e+07
#> TotFeASC          7.688476e+02   7.688476e+02   3.638801e+07   3.638801e+07
#> TotH2SoxidFe      1.909261e+00   1.909261e+00   7.064594e-07   7.105837e-07
#> TotH2SoxidMn      9.298853e-01   9.298853e-01   1.013465e-06   1.019381e-06
#> TotAgeMnox        2.365094e-01   2.365094e-01   4.955916e-01   4.950376e-01
#> TotAgeFeox        7.286289e-01   7.286289e-01   4.993288e-01   4.992824e-01
#> TotMnCO3prec      7.167273e-01   7.167273e-01   7.136831e-04   7.235036e-04
#> TotalMn           8.873666e-01   8.873666e-01   1.128465e-03   1.142716e-03
#> TotalMnO2         3.039111e+02   3.039111e+02   3.625077e+07   3.625077e+07
#> PartMnred         2.118412e-04   2.118412e-04   2.347379e-04   2.100068e-04
#> TotFeOxidMnASC    6.137837e-05   6.137837e-05   5.235894e-08   5.282645e-08
#> Cflux             4.566210e+02   4.566210e+02   1.000000e+01   1.000000e+01
#> FeOH3flux         1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00
#> CaPflux           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> w                 2.739726e-06   2.739726e-06   2.739726e-06   2.739726e-06
#> biotfac           1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00
#> irrfac            1.000000e+00   0.000000e+00   1.000000e+00   0.000000e+00
#> rFast             6.849315e-02   6.849315e-02   6.849315e-02   6.849315e-02
#> rSlow             1.369863e-04   1.369863e-04   1.369863e-04   1.369863e-04
#> pFast             9.000000e-01   9.000000e-01   9.000000e-01   9.000000e-01
#> MPBprod           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> gasflux           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> bwO2              3.000000e+02   3.000000e+02   3.000000e+02   3.000000e+02
#> bwNO3             1.000000e+01   1.000000e+01   1.000000e+01   1.000000e+01
#> bwNO2             0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> bwNH3             1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00
#> bwCH4             0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> bwFe              0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> bwH2S             0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> bwSO4             3.100000e+04   3.100000e+04   3.100000e+04   3.100000e+04
#> bwPO4             5.000000e-01   5.000000e-01   5.000000e-01   5.000000e-01
#> bwDIC             2.100000e+03   2.100000e+03   2.100000e+03   2.100000e+03
#> bwALK             2.400000e+03   2.400000e+03   2.400000e+03   2.400000e+03
#> Hwater            0.000000e+00   1.000000e+01   0.000000e+00   1.000000e+01
#> Ratefactor        1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00
#> bwMn              0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
#> MnO2flux          1.000000e+00   1.000000e+00   1.000000e+00   1.000000e+00

# 1-D outputs
  head(FESDIA1D(out))
#>            x       por     FDET     SDET       O2      NO3        NO2       NH3
#> 1 0.00498923 0.8934027 3405.716 3979.300 299.9200 10.01177 0.01168032 0.9851403
#> 2 0.01530360 0.8801069 3114.392 3946.070 299.7548 10.03681 0.03448383 0.9551439
#> 3 0.02631242 0.8664113 2845.251 3914.409 299.5789 10.06428 0.05732484 0.9239286
#> 4 0.03806243 0.8523376 2595.419 3883.986 299.3919 10.09439 0.08003927 0.8914989
#> 5 0.05060354 0.8379122 2362.751 3854.538 299.1935 10.12736 0.10244021 0.8578692
#> 6 0.06398901 0.8231666 2145.616 3825.852 298.9834 10.16340 0.12431765 0.8230660
#>        DIC           Fe    FeOH3          H2S   SO4          CH4       PO4
#> 1 2100.098 1.035468e-05 696.0885 1.913274e-09 31000 5.707679e-12 0.5011495
#> 2 2100.300 2.891200e-05 679.4919 5.633433e-09 31000 1.084172e-11 0.5035355
#> 3 2100.517 4.577695e-05 663.7130 9.313691e-09 31000 1.322379e-11 0.5060879
#> 4 2100.748 6.077689e-05 648.5875 1.289176e-08 31000 1.436682e-11 0.5088126
#> 5 2100.994 7.376661e-05 633.9870 1.629818e-08 31000 1.494457e-11 0.5117148
#> 6 2101.256 8.463450e-05 619.8085 1.945796e-08 31000 1.523563e-11 0.5147986
#>        FeP CaP Pads      ALK   FeOH3B           Mn     MnO2    MnO2B       TOC
#> 1 32.45068   0    0 2399.933 729152.0 6.003606e-05 467.4337 726629.3 0.5035448
#> 2 32.45109   0    0 2399.796 729142.9 1.813232e-04 450.8765 726620.1 0.5033890
#> 3 32.45159   0    0 2399.652 729135.0 3.072068e-04 435.1807 726612.0 0.5032446
#> 4 32.45221   0    0 2399.503 729128.1 4.371533e-04 420.1850 726605.1 0.5031101
#> 5 32.45294   0    0 2399.346 729122.3 5.704867e-04 405.7648 726599.0 0.5029843
#> 6 32.45381   0    0 2399.183 729117.4 7.063744e-04 391.8213 726593.9 0.5028663
#>   DICprodMin DINprodMin DIPprodMin O2prod  Oxicmin Denitrific    Feredmin
#> 1   27.89769   4.210972  0.2631858      0 27.86047 0.02339823 0.003134806
#> 2   29.13254   4.397365  0.2748353      0 29.09400 0.02449356 0.003216710
#> 3   30.13052   4.548003  0.2842502      0 30.09099 0.02539982 0.003270091
#> 4   30.88950   4.662566  0.2914104      0 30.84928 0.02611460 0.003295747
#> 5   31.40737   4.740736  0.2962960      0 31.36678 0.02663528 0.003294462
#> 6   31.68268   4.782291  0.2988932      0 31.64202 0.02695949 0.003267093
#>         BSRmin      Methmin   nitri1    nitri2     Anammox       Feoxid
#> 1 6.505098e-06 2.137764e-07 19.63733 0.2328302 0.001150676 0.0009316728
#> 2 6.826866e-06 2.243506e-07 19.03936 0.6873834 0.003293701 0.0025999537
#> 3 7.092388e-06 2.330764e-07 18.41710 1.1426825 0.005296406 0.0041141426
#> 4 7.300372e-06 2.399114e-07 17.77062 1.5954563 0.007135492 0.0054588336
#> 5 7.449602e-06 2.448155e-07 17.10023 2.0419793 0.008788030 0.0066211470
#> 6 7.539099e-06 2.477566e-07 16.40645 2.4780647 0.010232163 0.0075912934
#>        H2Soxid      CH4oxid          AOM      FeSprod    FePadsorp    FePdesorp
#> 1 2.869145e-10 4.621988e-08 5.308142e-12 1.981134e-17 0.0003488444 0.0005845613
#> 2 8.443244e-10 8.774617e-08 1.008280e-11 1.628738e-16 0.0003421483 0.0006144929
#> 3 1.395093e-09 1.069623e-07 1.229812e-11 4.263523e-16 0.0003358971 0.0006395516
#> 4 1.929845e-09 1.161353e-07 1.336114e-11 7.835214e-16 0.0003300095 0.0006596135
#> 5 2.438155e-09 1.207256e-07 1.389845e-11 1.202262e-15 0.0003244206 0.0006745563
#> 6 2.908804e-09 1.229904e-07 1.416913e-11 1.646815e-15 0.0003190766 0.0006842733
#>   CaPprod CaPdiss Padsorb H2Soxsurf CH4oxsurf O2distConsump   ALKprod
#> 1       0       0       0         0         0             0 -35.28530
#> 2       0       0       0         0         0             0 -33.91671
#> 3       0       0       0         0         0             0 -32.53312
#> 4       0       0       0         0         0             0 -31.13497
#> 5       0       0       0         0         0             0 -29.72311
#> 6       0       0       0         0         0             0 -28.29878
#>      DICprodCH4 MPBCprod MPBuptakeNO3 MPBuptakeNH3 MPBuptakePO4 MPBuptakeDIC
#> 1 -6.066301e-08        0            0            0            0            0
#> 2 -2.441907e-08        0            0            0            0            0
#> 3 -9.563583e-09        0            0            0            0            0
#> 4 -3.806986e-09        0            0            0            0            0
#> 5 -1.668301e-09        0            0            0            0            0
#> 6 -8.737746e-10        0            0            0            0            0
#>     MnredMin  sumMnO2 sumFeOH3  FeOxidMnASC H2SOxidFeOH3ASC H2SoxidMnO2ASC
#> 1 0.01068410 727096.7 729848.1 4.878698e-05    1.675679e-07   2.403881e-07
#> 2 0.01081958 727070.9 729822.4 1.362166e-04    4.933687e-07   7.077725e-07
#> 3 0.01085561 727047.2 729798.7 2.156674e-04    8.156544e-07   1.170114e-06
#> 4 0.01079771 727025.2 729776.7 2.863275e-04    1.128973e-06   1.619592e-06
#> 5 0.01065128 727004.8 729756.3 3.475138e-04    1.427244e-06   2.047483e-06
#> 6 0.01042188 726985.7 729737.2 3.987019e-04    1.703904e-06   2.444370e-06
#>        MnSprod       MnOxid   AgeFeOx  AgeMnOx    MnCO3prec
#> 1 1.148654e-18 1.555720e-05 1.0825568 2.180859 3.782448e-05
#> 2 1.021472e-17 4.696057e-05 1.0567458 2.103610 1.142500e-04
#> 3 2.861229e-17 7.951623e-05 1.0322064 2.030379 1.935879e-04
#> 4 5.635677e-17 1.130805e-04 1.0086833 1.960415 2.755047e-04
#> 5 9.297898e-17 1.474726e-04 0.9859766 1.893136 3.595768e-04
#> 6 1.374461e-16 1.824718e-04 0.9639262 1.828081 4.452821e-04
  head(FESDIA1D(defdyn, which = c("O2", "TOC")))
#>            x       por       O2       TOC
#> 1 0.00498923 0.8934027 297.6528 0.6618634
#> 2 0.01530360 0.8801069 292.7710 0.6547498
#> 3 0.02631242 0.8664113 287.5379 0.6481569
#> 4 0.03806243 0.8523376 281.9388 0.6420144
#> 5 0.05060354 0.8379122 275.9607 0.6362693
#> 6 0.06398901 0.8231666 269.5928 0.6308815

# and the pH
  mf <- par (mfrow = c(1,2))
  plot(x = FESDIApH(out), y = FESDIAdepth(dyna), ylim = c(10,0), type = "l")

  image2D(FESDIApH(dyna), y = FESDIAdepth(dyna), ylim = c(10,0))



  par(mfrow = mf)
```
