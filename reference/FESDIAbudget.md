# Budgets for the CNPDIA or MPBDIA model.

`FESDIAbudgetO2, FESDIAbudgetC, FESDIAbudgetN, FESDIAbudgetP` estimate
mass budgets from FESDIA model solutions.

## Usage

``` r
FESDIAbudgetO2(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetC(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetN(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetP(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetS(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetFe(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"))
```

## Arguments

- out :

  an output object returned by
  [FESDIAsolve](https://stanleesocca.github.io/FESDIA/reference/FESDIAsolve.md)
  or
  [FESDIAdyna](https://stanleesocca.github.io/FESDIA/reference/FESDIAdyna.md).

- which :

  if not `NULL`, a vector with names of the items to return.

- ... :

  unused.

## Value

`FESDIAbudgetx` returns the element budget (C, N, P, S, Fe, O2) of the
solution as a `list`, with the following items.

- `Fluxes`, the boundary fluxes at the surface and bottom of the
  sediment, the perturbation fluxes (only when the model was solved with
  [FESDIAperturb](https://stanleesocca.github.io/FESDIA/reference/FESDIAperturb.md))
  and the net input. Positive fluxes are directed into the sediment for
  the surface, and out of the sediment at the bottom. Negative
  perturbation fluxes are directed out of the sediment. For dynamic
  runs, fluxes are averaged over the simulation period.

- `Rates`, the integrated process rates, in nmol/cm2/d

- `dC`, the rate of change of integrated values of state vairables, in
  nmol/cm2/d, defined as (\[concentration at the end\] - \[concentration
  at beginning\])/\[length of simulation\]

- `Losses`, the total amount lost from the system (burial, and removal
  e.g. N2 production for N-budget)

- `Delta`, the difference between total fluxes in and fluxes out, i.e.
  the deviation from steady-state.

- `Fluxmat`, the flux matrix (rows: from, columns:to). For N, C, Fe, S
  and P budget the column sums - row sums is equal to the rate of
  change. For the O2 budget, this only applies to O2.

If more than one FESDIA object is passed to these functions, a matrix is
returned, one column for each object.

For dynamic runs, the budget is taken over the mean of the run; Delta is
then the integrated mean rate of change; for steady-state runs, Delta
should be very small.

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
# some runs to work with  
  defsteady <- FESDIAsolve()
  defdyn    <- FESDIAdyna()

# altered steady-state run
  out <- FESDIAsolve(parms = list(Cflux = 1000))
  cbind(default = FESDIAparms(), altered = FESDIAparms(out))
#>               default.parms default.units
#> Cflux          4.566210e+02   nmolC/cm2/d
#> pFast          9.000000e-01             -
#> FeOH3flux      1.000000e+00    nmol/cm2/d
#> CaPflux        0.000000e+00    nmol/cm2/d
#> rFast          6.849315e-02            /d
#> rSlow          1.369863e-04            /d
#> NCrFdet        1.509434e-01     molN/molC
#> NCrSdet        1.509434e-01     molN/molC
#> PCrFdet        9.433962e-03     molP/molC
#> PCrSdet        9.433962e-03     molP/molC
#> BCupLiq        2.000000e+00             -
#> BCdownLiq      3.000000e+00             -
#> O2bw           3.000000e+02       mmol/m3
#> NO3bw          1.000000e+01       mmol/m3
#> NO2bw          0.000000e+00       mmol/m3
#> NH3bw          1.000000e+00       mmol/m3
#> CH4bw          0.000000e+00       mmol/m3
#> PO4bw          5.000000e-01       mmol/m3
#> DICbw          2.100000e+03       mmol/m3
#> Febw           0.000000e+00       mmol/m3
#> H2Sbw          0.000000e+00       mmol/m3
#> SO4bw          3.100000e+04       mmol/m3
#> ALKbw          2.400000e+03       mmol/m3
#> Mnbw           0.000000e+00       mmol/m3
#> O2dw                     NA       mmol/m3
#> NO3dw                    NA       mmol/m3
#> NO2dw                    NA       mmol/m3
#> NH3dw                    NA       mmol/m3
#> CH4dw                    NA       mmol/m3
#> PO4dw                    NA       mmol/m3
#> DICdw                    NA       mmol/m3
#> Fedw                     NA       mmol/m3
#> H2Sdw                    NA       mmol/m3
#> SO4dw                    NA       mmol/m3
#> ALKdw                    NA       mmol/m3
#> Mndw                     NA       mmol/m3
#> w              2.739726e-06          cm/d
#> biot           2.739726e-03         cm2/d
#> biotdepth      5.000000e+00            cm
#> biotatt        1.000000e+00           /cm
#> irr            0.000000e+00            /d
#> irrdepth       5.000000e+00            cm
#> irratt         1.000000e+00            cm
#> gasflux        0.000000e+00          cm/d
#> NH3Ads         1.300000e+00             -
#> rnitri1        2.000000e+01            /d
#> rnitri2        2.000000e+01            /d
#> ranammox       1.000000e-01  /(mmol/m3)/d
#> ksO2nitri      1.000000e+00     mmolO2/m3
#> ksO2oxic       3.000000e+00     mmolO2/m3
#> ksNO3denit     3.000000e+01    mmolNO3/m3
#> kinO2denit     1.000000e+00     mmolO2/m3
#> kinNO3anox     1.000000e+00    mmolNO3/m3
#> kinO2anox      1.000000e-03     mmolO2/m3
#> temperature    1.000000e+01           dgC
#> salinity       3.500000e+01           psu
#> TOC0           5.000000e-01             %
#> rFePadsorp     1.000000e-06            /d
#> rCaPprod       0.000000e+00            /d
#> rCaPdiss       0.000000e+00            /d
#> CPrCaP         2.869565e-01       mol/mol
#> rPads          0.000000e+00            /d
#> rPdes          0.000000e+00            /d
#> maxPads        1.000000e+03 mmolP/m3solid
#> ksFeOH3        1.250000e+04  mmolFeOH3/m3
#> kinFeOH3       1.250000e+04  mmolFeOH3/m3
#> ksSO4BSR       1.600000e+03      mmolS/m3
#> kinSO4Met      1.000000e+03      mmolS/m3
#> rFeox          3.000000e-01  /(mmol/m3)/d
#> rH2Sox         5.000000e-04  /(mmol/m3)/d
#> rFeS           1.000000e-03  /(mmol/m3)/d
#> rCH4ox         2.700000e+01  /(mmol/m3)/d
#> rAOM           3.000000e-05  /(mmol/m3)/d
#> rSurfH2Sox     0.000000e+00            /d
#> rSurfCH4ox     0.000000e+00            /d
#> ksSurfALK      3.000000e+03       mmol/m3
#> ksO2reox       1.000000e+00     mmolO2/m3
#> ODUoxdepth     5.000000e+00            cm
#> ODUoxatt       1.000000e+00           /cm
#> por0           9.000000e-01             -
#> pordeep        5.000000e-01             -
#> porcoeff       3.000000e-01            cm
#> formationtype  1.000000e+00             -
#> dilution       0.000000e+00            /d
#> Hwater         1.000000e+01            cm
#> Cfall          1.000000e+02          cm/d
#> FePfall        1.000000e+02          cm/d
#> FeOH3fall      1.000000e+02          cm/d
#> CaPfall        1.000000e+02          cm/d
#> addalk         1.000000e+00             -
#> MPBprod        0.000000e+00     mmol/m3/d
#> kMPB           4.000000e+00           /cm
#> kDINupt        1.000000e-02       mmol/m3
#> kPO4upt        1.000000e-03       mmol/m3
#> kDICupt        1.000000e+00       mmol/m3
#> rH2Sfeox       1.200000e-04    cm3/nmol/d
#> MnO2flux       1.000000e+00    nmol/cm2/d
#> rAgeFeox       1.555200e-03    cm3/nmol/d
#> rMnOxid        8.640000e-04    cm3/nmol/d
#> rH2SMnox       1.728000e-04    cm3/nmol/d
#> rAgeMnox       4.665600e-03    cm3/nmol/d
#> rMnFe          6.480000e-06    cm3/nmol/d
#> rMnS           1.000000e-05    cm3/nmol/d
#> rMnCO3prec     3.000000e-04    cm3/nmol/d
#> rFeCO3prec     3.000000e-04    cm3/nmol/d
#> ksMnO2         2.600000e+03       mmol/m3
#> pFastFeOx      5.000000e-01             -
#> pFastMnOx      5.000000e-01             -
#> kinMnO2        2.600000e+03       mmol/m3
#> isDICcorr      0.000000e+00             -
#>                                                    default.description
#> Cflux                                       total organic C deposition
#> pFast                                         part FDET in carbon flux
#> FeOH3flux                                     deposition rate of FeOH3
#> CaPflux                                         deposition rate of CaP
#> rFast                                                  decay rate FDET
#> rSlow                                                  decay rate SDET
#> NCrFdet                                                  NC ratio FDET
#> NCrSdet                                                  NC ratio SDET
#> PCrFdet                                                  PC ratio FDET
#> PCrSdet                                                  PC ratio SDET
#> BCupLiq                   upper boundary liq. 1:flux, 2:conc, 3:0-grad
#> BCdownLiq                 lower boundary liq. 1:flux, 2:conc, 3:0-grad
#> O2bw                         upper boundary O2  -if BC=1: flux, 2:conc
#> NO3bw                        upper boundary NO3 -if BC=1: flux, 2:conc
#> NO2bw                        upper boundary NO2 -if BC=1: flux, 2:conc
#> NH3bw                        upper boundary NH3 -if BC=1: flux, 2:conc
#> CH4bw                        upper boundary CH4 -if BC=1: flux, 2:conc
#> PO4bw                        upper boundary PO4 -if BC=1: flux, 2:conc
#> DICbw                        upper boundary DIC -if BC=1: flux, 2:conc
#> Febw                        upper boundary Fe2+ -if BC=1: flux, 2:conc
#> H2Sbw                        upper boundary H2S -if BC=1: flux, 2:conc
#> SO4bw                        upper boundary SO4 -if BC=1: flux, 2:conc
#> ALKbw                 upper boundary alkalinity -if BC=1: flux, 2:conc
#> Mnbw                   upper boundary Manganese -if BC=1: flux, 2:conc
#> O2dw                         lower boundary O2  -if BC=1: flux, 2:conc
#> NO3dw                        lower boundary NO3 -if BC=1: flux, 2:conc
#> NO2dw                        lower boundary NO2 -if BC=1: flux, 2:conc
#> NH3dw                        lower boundary NH3 -if BC=1: flux, 2:conc
#> CH4dw                        lower boundary CH3 -if BC=1: flux, 2:conc
#> PO4dw                        lower boundary PO4 -if BC=1: flux, 2:conc
#> DICdw                        lower boundary DIC -if BC=1: flux, 2:conc
#> Fedw                        lower boundary Fe2+ -if BC=1: flux, 2:conc
#> H2Sdw                        lower boundary H2S -if BC=1: flux, 2:conc
#> SO4dw                        lower boundary SO4 -if BC=1: flux, 2:conc
#> ALKdw                 lower boundary alkalinity -if BC=1: flux, 2:conc
#> Mndw                   lower boundary Mangenese -if BC=1: flux, 2:conc
#> w                                                       advection rate
#> biot                                          bioturbation coefficient
#> biotdepth                                         depth of mixed layer
#> biotatt                              attenuation coeff below biotdepth
#> irr                                                bio-irrigation rate
#> irrdepth                                      depth of irrigated layer
#> irratt                                attenuation coeff below irrdepth
#> gasflux                                  piston velocity for dry flats
#> NH3Ads                                       Adsorption coeff ammonium
#> rnitri1                           Max nitrification rate step1 (NH3ox)
#> rnitri2                           Max nitrification rate step2 (NO2ox)
#> ranammox                                                  Anammox rate
#> ksO2nitri                                 half-sat O2 in nitrification
#> ksO2oxic                            half-sat O2 in oxic mineralisation
#> ksNO3denit                             half-sat NO3 in denitrification
#> kinO2denit                           half-sat O2 inhib denitrification
#> kinNO3anox                              half-sat NO3 inhib anoxic degr
#> kinO2anox                                 half-sat O2 inhib anoxic min
#> temperature                                                temperature
#> salinity                                                      salinity
#> TOC0                                            refractory Carbon conc
#> rFePadsorp                                         rate FeP adsorption
#> rCaPprod                                           rate CaP production
#> rCaPdiss                                          rate CaP dissolution
#> CPrCaP                                                 C:Pratio in CaP
#> rPads                                              adsorption rate PO4
#> rPdes                                    desorption rate of adsorbed P
#> maxPads                                   Max adsorbed P concentration
#> ksFeOH3                          half-sat FeOH3 conc in iron reduction
#> kinFeOH3                         half-sat FeOH3 inhibition S reduction
#> ksSO4BSR                       half-sat SO4 conc in sulphate reduction
#> kinSO4Met                       half-sat SO4 inhibition methanogenesis
#> rFeox                                            Max rate Fe oxidation
#> rH2Sox                                          Max rate H2S oxidation
#> rFeS                                       maximum rate FeS production
#> rCH4ox                                  Max rate CH4 oxidation with O2
#> rAOM                              Max rate anaerobic oxidation Methane
#> rSurfH2Sox                           Max rate H2S oxidation with BW O2
#> rSurfCH4ox                           Max rate CH4 oxidation with BW O2
#> ksSurfALK        half-sat Alkalinity in oxidation of H2S/CH4 with bwO2
#> ksO2reox             half-sat Oxygen in oxidation of H2S/CH4 with bwO2
#> ODUoxdepth                      Max depth H2S/CH4 oxidation with BW O2
#> ODUoxatt                               depth attenuation ODU oxidation
#> por0                                                  surface porosity
#> pordeep                                                  deep porosity
#> porcoeff                                    porosity decay coefficient
#> formationtype            formationfactor, 1=sand,2=fine sand,3=general
#> dilution                           relaxation towards background conc 
#> Hwater                                       height of water over core
#> Cfall                             fall speed of organic C (FDET, SDET)
#> FePfall                                              fall speed of FeP
#> FeOH3fall                                          fall speed of FeOH3
#> CaPfall                                              fall speed of CaP
#> addalk                                            solve for alkalinity
#> MPBprod                                    maximal MPB production rate
#> kMPB                          sedimentary light extinction coefficient
#> kDINupt                                    DIN limitation constant MPB
#> kPO4upt                                      P limitation constant MPB
#> kDICupt                                      C limitation constant MPB
#> rH2Sfeox      Rate of Sulphide-mediated iron reduction (oxyhydr)oxides
#> MnO2flux                                             Flux of Mn Oxides
#> rAgeFeox                                              Ageing of FeOx=A
#> rMnOxid                            Rate of Reoxidation of Mn by Oxygen
#> rH2SMnox                            Rate of Reoxidation of H2S by MnOx
#> rAgeMnox                                              Ageing of MnOx=A
#> rMnFe                              Rate of Reoxidation of Fe with MnOx
#> rMnS                                             Rate of MnS formation
#> rMnCO3prec                                     Rate of MnCO3 formation
#> rFeCO3prec                                     Rate of FeCO3 formation
#> ksMnO2                             Half saturation constant for Mn red
#> pFastFeOx                                  fraction fast FeOH3 in flux
#> pFastMnOx                                   fraction fast MnO2 in flux
#> kinMnO2                           Inhibition of iron and other by MnO2
#> isDICcorr      DIC correction for Calcite precip - Rassmann et al 2020
#>               altered.parms altered.units
#> Cflux          1.000000e+03   nmolC/cm2/d
#> pFast          9.000000e-01             -
#> FeOH3flux      1.000000e+00    nmol/cm2/d
#> CaPflux        0.000000e+00    nmol/cm2/d
#> rFast          6.849315e-02            /d
#> rSlow          1.369863e-04            /d
#> NCrFdet        1.509434e-01     molN/molC
#> NCrSdet        1.509434e-01     molN/molC
#> PCrFdet        9.433962e-03     molP/molC
#> PCrSdet        9.433962e-03     molP/molC
#> BCupLiq        2.000000e+00             -
#> BCdownLiq      3.000000e+00             -
#> O2bw           3.000000e+02       mmol/m3
#> NO3bw          1.000000e+01       mmol/m3
#> NO2bw          0.000000e+00       mmol/m3
#> NH3bw          1.000000e+00       mmol/m3
#> CH4bw          0.000000e+00       mmol/m3
#> PO4bw          5.000000e-01       mmol/m3
#> DICbw          2.100000e+03       mmol/m3
#> Febw           0.000000e+00       mmol/m3
#> H2Sbw          0.000000e+00       mmol/m3
#> SO4bw          3.100000e+04       mmol/m3
#> ALKbw          2.400000e+03       mmol/m3
#> Mnbw           0.000000e+00       mmol/m3
#> O2dw                     NA       mmol/m3
#> NO3dw                    NA       mmol/m3
#> NO2dw                    NA       mmol/m3
#> NH3dw                    NA       mmol/m3
#> CH4dw                    NA       mmol/m3
#> PO4dw                    NA       mmol/m3
#> DICdw                    NA       mmol/m3
#> Fedw                     NA       mmol/m3
#> H2Sdw                    NA       mmol/m3
#> SO4dw                    NA       mmol/m3
#> ALKdw                    NA       mmol/m3
#> Mndw                     NA       mmol/m3
#> w              2.739726e-06          cm/d
#> biot           2.739726e-03         cm2/d
#> biotdepth      5.000000e+00            cm
#> biotatt        1.000000e+00           /cm
#> irr            0.000000e+00            /d
#> irrdepth       5.000000e+00            cm
#> irratt         1.000000e+00            cm
#> gasflux        0.000000e+00          cm/d
#> NH3Ads         1.300000e+00             -
#> rnitri1        2.000000e+01            /d
#> rnitri2        2.000000e+01            /d
#> ranammox       1.000000e-01  /(mmol/m3)/d
#> ksO2nitri      1.000000e+00     mmolO2/m3
#> ksO2oxic       3.000000e+00     mmolO2/m3
#> ksNO3denit     3.000000e+01    mmolNO3/m3
#> kinO2denit     1.000000e+00     mmolO2/m3
#> kinNO3anox     1.000000e+00    mmolNO3/m3
#> kinO2anox      1.000000e-03     mmolO2/m3
#> temperature    1.000000e+01           dgC
#> salinity       3.500000e+01           psu
#> TOC0           5.000000e-01             %
#> rFePadsorp     1.000000e-06            /d
#> rCaPprod       0.000000e+00            /d
#> rCaPdiss       0.000000e+00            /d
#> CPrCaP         2.869565e-01       mol/mol
#> rPads          0.000000e+00            /d
#> rPdes          0.000000e+00            /d
#> maxPads        1.000000e+03 mmolP/m3solid
#> ksFeOH3        1.250000e+04  mmolFeOH3/m3
#> kinFeOH3       1.250000e+04  mmolFeOH3/m3
#> ksSO4BSR       1.600000e+03      mmolS/m3
#> kinSO4Met      1.000000e+03      mmolS/m3
#> rFeox          3.000000e-01  /(mmol/m3)/d
#> rH2Sox         5.000000e-04  /(mmol/m3)/d
#> rFeS           1.000000e-03  /(mmol/m3)/d
#> rCH4ox         2.700000e+01  /(mmol/m3)/d
#> rAOM           3.000000e-05  /(mmol/m3)/d
#> rSurfH2Sox     0.000000e+00            /d
#> rSurfCH4ox     0.000000e+00            /d
#> ksSurfALK      3.000000e+03       mmol/m3
#> ksO2reox       1.000000e+00     mmolO2/m3
#> ODUoxdepth     5.000000e+00            cm
#> ODUoxatt       1.000000e+00           /cm
#> por0           9.000000e-01             -
#> pordeep        5.000000e-01             -
#> porcoeff       3.000000e-01            cm
#> formationtype  1.000000e+00             -
#> dilution       0.000000e+00            /d
#> Hwater         1.000000e+01            cm
#> Cfall          1.000000e+02          cm/d
#> FePfall        1.000000e+02          cm/d
#> FeOH3fall      1.000000e+02          cm/d
#> CaPfall        1.000000e+02          cm/d
#> addalk         1.000000e+00             -
#> MPBprod        0.000000e+00     mmol/m3/d
#> kMPB           4.000000e+00           /cm
#> kDINupt        1.000000e-02       mmol/m3
#> kPO4upt        1.000000e-03       mmol/m3
#> kDICupt        1.000000e+00       mmol/m3
#> rH2Sfeox       1.200000e-04    cm3/nmol/d
#> MnO2flux       1.000000e+00    nmol/cm2/d
#> rAgeFeox       1.555200e-03    cm3/nmol/d
#> rMnOxid        8.640000e-04    cm3/nmol/d
#> rH2SMnox       1.728000e-04    cm3/nmol/d
#> rAgeMnox       4.665600e-03    cm3/nmol/d
#> rMnFe          6.480000e-06    cm3/nmol/d
#> rMnS           1.000000e-05    cm3/nmol/d
#> rMnCO3prec     3.000000e-04    cm3/nmol/d
#> rFeCO3prec     3.000000e-04    cm3/nmol/d
#> ksMnO2         2.600000e+03       mmol/m3
#> pFastFeOx      5.000000e-01             -
#> pFastMnOx      5.000000e-01             -
#> kinMnO2        2.600000e+03       mmol/m3
#> isDICcorr      0.000000e+00             -
#>                                                    altered.description
#> Cflux                                       total organic C deposition
#> pFast                                         part FDET in carbon flux
#> FeOH3flux                                     deposition rate of FeOH3
#> CaPflux                                         deposition rate of CaP
#> rFast                                                  decay rate FDET
#> rSlow                                                  decay rate SDET
#> NCrFdet                                                  NC ratio FDET
#> NCrSdet                                                  NC ratio SDET
#> PCrFdet                                                  PC ratio FDET
#> PCrSdet                                                  PC ratio SDET
#> BCupLiq                   upper boundary liq. 1:flux, 2:conc, 3:0-grad
#> BCdownLiq                 lower boundary liq. 1:flux, 2:conc, 3:0-grad
#> O2bw                         upper boundary O2  -if BC=1: flux, 2:conc
#> NO3bw                        upper boundary NO3 -if BC=1: flux, 2:conc
#> NO2bw                        upper boundary NO2 -if BC=1: flux, 2:conc
#> NH3bw                        upper boundary NH3 -if BC=1: flux, 2:conc
#> CH4bw                        upper boundary CH4 -if BC=1: flux, 2:conc
#> PO4bw                        upper boundary PO4 -if BC=1: flux, 2:conc
#> DICbw                        upper boundary DIC -if BC=1: flux, 2:conc
#> Febw                        upper boundary Fe2+ -if BC=1: flux, 2:conc
#> H2Sbw                        upper boundary H2S -if BC=1: flux, 2:conc
#> SO4bw                        upper boundary SO4 -if BC=1: flux, 2:conc
#> ALKbw                 upper boundary alkalinity -if BC=1: flux, 2:conc
#> Mnbw                   upper boundary Manganese -if BC=1: flux, 2:conc
#> O2dw                         lower boundary O2  -if BC=1: flux, 2:conc
#> NO3dw                        lower boundary NO3 -if BC=1: flux, 2:conc
#> NO2dw                        lower boundary NO2 -if BC=1: flux, 2:conc
#> NH3dw                        lower boundary NH3 -if BC=1: flux, 2:conc
#> CH4dw                        lower boundary CH3 -if BC=1: flux, 2:conc
#> PO4dw                        lower boundary PO4 -if BC=1: flux, 2:conc
#> DICdw                        lower boundary DIC -if BC=1: flux, 2:conc
#> Fedw                        lower boundary Fe2+ -if BC=1: flux, 2:conc
#> H2Sdw                        lower boundary H2S -if BC=1: flux, 2:conc
#> SO4dw                        lower boundary SO4 -if BC=1: flux, 2:conc
#> ALKdw                 lower boundary alkalinity -if BC=1: flux, 2:conc
#> Mndw                   lower boundary Mangenese -if BC=1: flux, 2:conc
#> w                                                       advection rate
#> biot                                          bioturbation coefficient
#> biotdepth                                         depth of mixed layer
#> biotatt                              attenuation coeff below biotdepth
#> irr                                                bio-irrigation rate
#> irrdepth                                      depth of irrigated layer
#> irratt                                attenuation coeff below irrdepth
#> gasflux                                  piston velocity for dry flats
#> NH3Ads                                       Adsorption coeff ammonium
#> rnitri1                           Max nitrification rate step1 (NH3ox)
#> rnitri2                           Max nitrification rate step2 (NO2ox)
#> ranammox                                                  Anammox rate
#> ksO2nitri                                 half-sat O2 in nitrification
#> ksO2oxic                            half-sat O2 in oxic mineralisation
#> ksNO3denit                             half-sat NO3 in denitrification
#> kinO2denit                           half-sat O2 inhib denitrification
#> kinNO3anox                              half-sat NO3 inhib anoxic degr
#> kinO2anox                                 half-sat O2 inhib anoxic min
#> temperature                                                temperature
#> salinity                                                      salinity
#> TOC0                                            refractory Carbon conc
#> rFePadsorp                                         rate FeP adsorption
#> rCaPprod                                           rate CaP production
#> rCaPdiss                                          rate CaP dissolution
#> CPrCaP                                                 C:Pratio in CaP
#> rPads                                              adsorption rate PO4
#> rPdes                                    desorption rate of adsorbed P
#> maxPads                                   Max adsorbed P concentration
#> ksFeOH3                          half-sat FeOH3 conc in iron reduction
#> kinFeOH3                         half-sat FeOH3 inhibition S reduction
#> ksSO4BSR                       half-sat SO4 conc in sulphate reduction
#> kinSO4Met                       half-sat SO4 inhibition methanogenesis
#> rFeox                                            Max rate Fe oxidation
#> rH2Sox                                          Max rate H2S oxidation
#> rFeS                                       maximum rate FeS production
#> rCH4ox                                  Max rate CH4 oxidation with O2
#> rAOM                              Max rate anaerobic oxidation Methane
#> rSurfH2Sox                           Max rate H2S oxidation with BW O2
#> rSurfCH4ox                           Max rate CH4 oxidation with BW O2
#> ksSurfALK        half-sat Alkalinity in oxidation of H2S/CH4 with bwO2
#> ksO2reox             half-sat Oxygen in oxidation of H2S/CH4 with bwO2
#> ODUoxdepth                      Max depth H2S/CH4 oxidation with BW O2
#> ODUoxatt                               depth attenuation ODU oxidation
#> por0                                                  surface porosity
#> pordeep                                                  deep porosity
#> porcoeff                                    porosity decay coefficient
#> formationtype            formationfactor, 1=sand,2=fine sand,3=general
#> dilution                           relaxation towards background conc 
#> Hwater                                       height of water over core
#> Cfall                             fall speed of organic C (FDET, SDET)
#> FePfall                                              fall speed of FeP
#> FeOH3fall                                          fall speed of FeOH3
#> CaPfall                                              fall speed of CaP
#> addalk                                            solve for alkalinity
#> MPBprod                                    maximal MPB production rate
#> kMPB                          sedimentary light extinction coefficient
#> kDINupt                                    DIN limitation constant MPB
#> kPO4upt                                      P limitation constant MPB
#> kDICupt                                      C limitation constant MPB
#> rH2Sfeox      Rate of Sulphide-mediated iron reduction (oxyhydr)oxides
#> MnO2flux                                             Flux of Mn Oxides
#> rAgeFeox                                              Ageing of FeOx=A
#> rMnOxid                            Rate of Reoxidation of Mn by Oxygen
#> rH2SMnox                            Rate of Reoxidation of H2S by MnOx
#> rAgeMnox                                              Ageing of MnOx=A
#> rMnFe                              Rate of Reoxidation of Fe with MnOx
#> rMnS                                             Rate of MnS formation
#> rMnCO3prec                                     Rate of MnCO3 formation
#> rFeCO3prec                                     Rate of FeCO3 formation
#> ksMnO2                             Half saturation constant for Mn red
#> pFastFeOx                                  fraction fast FeOH3 in flux
#> pFastMnOx                                   fraction fast MnO2 in flux
#> kinMnO2                           Inhibition of iron and other by MnO2
#> isDICcorr      DIC correction for Calcite precip - Rassmann et al 2020

# budgets
  B1 <- FESDIAbudgetO2(out)  
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
  B1$Fluxmat
#>        Ext       O2     NO2     NO3      DIC       SO4     FeOH3        MnO2
#> Ext      0 975.5627   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> O2       0   0.0000 117.195 22.3567 834.9715 0.5957938 0.4392728 0.004380027
#> NO2      0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> NO3      0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> DIC      0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> SO4      0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> FeOH3    0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> MnO2     0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#> Burial   0   0.0000   0.000  0.0000   0.0000 0.0000000 0.0000000 0.000000000
#>        Burial
#> Ext         0
#> O2          0
#> NO2         0
#> NO3         0
#> DIC         0
#> SO4         0
#> FeOH3       0
#> MnO2        0
#> Burial      0
  colSums(B1$Fluxmat) - rowSums(B1$Fluxmat)  # Small only for O2
#>           Ext            O2           NO2           NO3           DIC 
#> -9.755627e+02 -1.136868e-12  1.171950e+02  2.235670e+01  8.349715e+02 
#>           SO4         FeOH3          MnO2        Burial 
#>  5.957938e-01  4.392728e-01  4.380027e-03  0.000000e+00 
  FESDIAbudgetO2(out, defsteady)$Rates  
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#>                            [,1]         [,2]
#> Nitrification      1.395517e+02  78.36444216
#> FeOxidation        4.392728e-01   0.75792271
#> MnOxidation        4.380027e-03   0.01518014
#> H2Soxidation       5.957938e-01   0.73687837
#> CH4oxidation       1.230501e+00   0.03590063
#> H2Soxid.dist       0.000000e+00   0.00000000
#> CH4oxid.dist       0.000000e+00   0.00000000
#> OxicMineralisation 8.337410e+02 427.04955985
#> MPBO2production    0.000000e+00   0.00000000
#> MPBO2respiration   0.000000e+00   0.00000000
#> Total              9.755627e+02 506.95988386
  
  B2 <- FESDIAbudgetC(out)
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
  colSums(B2$Fluxmat) - rowSums(B2$Fluxmat)  # Small for all states
#>           Ext           DET           DIC           CaP           CH4 
#> -1.258111e+00  2.273737e-13  7.525396e-01  0.000000e+00  1.076012e-07 
#>           MPB         MnCO3         CaCO3          ARAG        Burial 
#>  0.000000e+00  4.990915e-01  0.000000e+00  0.000000e+00  6.479902e-03 
  
  FESDIAbudgetC(out,defsteady)
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
#> $Fluxes
#>                        [,1]           [,2]
#> FDETsurf       9.000000e+02   4.109589e+02
#> FDETdeep       0.000000e+00  1.001633e-191
#> FDETperturb    0.000000e+00   0.000000e+00
#> FDETnet        9.000000e+02   4.109589e+02
#> SDETsurf       1.000000e+02   4.566210e+01
#> SDETdeep       2.558897e-71   1.168446e-71
#> SDETperturb    0.000000e+00   0.000000e+00
#> SDETnet        1.000000e+02   4.566210e+01
#> DICsurf       -9.987419e+02  -4.550175e+02
#> DICdeep        6.479902e-03   4.502289e-03
#> DICperturb     0.000000e+00   0.000000e+00
#> DICnet        -9.987484e+02  -4.550220e+02
#> CH4surf       -1.076012e-07  -4.391394e-08
#> CH4deep        1.085043e-43   1.702274e-44
#> CH4perturb     0.000000e+00   0.000000e+00
#> CH4net        -1.076012e-07  -4.391394e-08
#> CinCaPsurf     0.000000e+00   0.000000e+00
#> CinCaPdeep     0.000000e+00  3.337892e-141
#> CinCaPperturb  0.000000e+00   0.000000e+00
#> CinCaPnet      0.000000e+00 -3.337892e-141
#> CaCO3surf      0.000000e+00   0.000000e+00
#> CaCO3deep      0.000000e+00   0.000000e+00
#> CaCO3perturb   0.000000e+00   0.000000e+00
#> CaCO3net       0.000000e+00   0.000000e+00
#> ARAGsurf       0.000000e+00   0.000000e+00
#> ARAGdeep       0.000000e+00   0.000000e+00
#> ARAGperturb    0.000000e+00   0.000000e+00
#> ARAGnet        0.000000e+00   0.000000e+00
#> Totalsurf      1.258111e+00   1.603499e+00
#> Totaldeep      6.479902e-03   4.502289e-03
#> Totalperturb   0.000000e+00   0.000000e+00
#> Totalnet       1.251631e+00   1.598997e+00
#> 
#> $Rates
#>                             [,1]         [,2]
#> OxicMineralisation   833.7410098 427.04955985
#> Denitrification       34.5749182   9.43748998
#> ManganeseReduction     0.2150795   0.09673114
#> IronReduction          0.4556142   0.12539842
#> SulphateReduction    126.7648882  19.27246892
#> Methanogenesis         4.2484901   0.63935626
#> TotalMineralisation 1000.0000000 456.62100457
#> CH4oxidation           0.6152507   0.01795032
#> MnCO3precitation       0.4990915   0.71672732
#> CH4oxid.dist           0.0000000   0.00000000
#> CH4oxidAOM             1.5089943   0.30172777
#> MPBDICuptake           0.0000000   0.00000000
#> MPBFDETproduction      0.0000000   0.00000000
#> MPBResp                0.0000000   0.00000000
#> CaPprecipitation       0.0000000   0.00000000
#> CaPdissolution         0.0000000   0.00000000
#> CaCO3dissolution       0.0000000   0.00000000
#> ARAGdissolution        0.0000000   0.00000000
#> CaCO3production        0.0000000   0.00000000
#> 
#> $Losses
#>             [,1]        [,2]
#> [1,] 0.006479902 0.004502289
#> 
#> $dC
#>                 [,1]           [,2]
#> Ext     1.258111e+00   1.603500e+00
#> DET    -2.273737e-13  -1.136868e-13
#> DIC    -7.525396e-01  -8.822699e-01
#> CaP     0.000000e+00  3.337892e-141
#> CH4    -1.076012e-07  -4.391394e-08
#> MPB     0.000000e+00   0.000000e+00
#> MnCO3  -4.990915e-01  -7.167273e-01
#> CaCO3   0.000000e+00   0.000000e+00
#> ARAG    0.000000e+00   0.000000e+00
#> Burial -6.479902e-03  -4.502289e-03
#> sum    -5.741935e-15   7.811460e-15
#> 
#> $Delta
#>          [,1]     [,2]
#> [1,] 1.251631 1.598997
#> 
  FESDIAbudgetC(out,defsteady,defdyn, which = "Rates")
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
#> Warning: number of columns of result, 7, is not a multiple of vector length 4 of arg 2
#> Warning: number of rows of result is not a multiple of vector length (arg 2)
#>                             [,1]         [,2]         [,3]
#> OxicMineralisation   833.7410098 427.04955985 427.04955985
#> Denitrification       34.5749182   9.43748998   9.43748998
#> ManganeseReduction     0.2150795   0.09673114   0.09673114
#> IronReduction          0.4556142   0.12539842   0.12539842
#> SulphateReduction    126.7648882  19.27246892  19.27246892
#> Methanogenesis         4.2484901   0.63935626   0.63935626
#> TotalMineralisation 1000.0000000 456.62100457 456.62100457
#> CH4oxidation           0.6152507   0.01795032   0.01795032
#> MnCO3precitation       0.4990915   0.71672732   0.71672732
#> CH4oxid.dist           0.0000000   0.00000000   0.00000000
#> CH4oxidAOM             1.5089943   0.30172777   0.30172777
#> MPBDICuptake           0.0000000   0.00000000   0.00000000
#> MPBFDETproduction      0.0000000   0.00000000   0.00000000
#> MPBResp                0.0000000   0.00000000   0.00000000
#> CaPprecipitation       0.0000000   0.00000000   0.00000000
#> CaPdissolution         0.0000000   0.00000000   0.00000000
#> CaCO3dissolution       0.0000000   0.00000000   0.00000000
#> ARAGdissolution        0.0000000   0.00000000   0.00000000
#> CaCO3production        0.0000000   0.00000000   0.00000000
```
