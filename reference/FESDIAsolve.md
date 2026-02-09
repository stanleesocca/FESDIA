# Steady-state solution of the FESDIA model, calculating C, N, P, O2, Fe and S dynamics in the sediment.

`FESDIAsolve` finds the steady-state solution of the FESDIA model; pH
can be calculated afterwards, ie ignoring carbonate dynamics.

## Usage

``` r
FESDIAsolve (parms = list(), yini = NULL, gridtype = 1, Grid = NULL, 
      porosity = NULL, bioturbation = NULL, irrigation = NULL,   
      surface = NULL, diffusionfactor = NULL, 
      dynamicbottomwater = FALSE, ratefactor = NULL, 
      calcpH = FALSE, verbose = FALSE, 
      method = NULL, times = c(0, 1e+06), ...)
```

## Arguments

- parms:

  A list with parameter values.Available parameters can be listed using
  function
  [FESDIAparms](https://stanleesocca.github.io/FESDIA/reference/FESDIAparms.md).

  See details of
  [FESDIAparms](https://stanleesocca.github.io/FESDIA/reference/FESDIAparms.md).

- gridtype :

  Type of grid: `1` for cartesian, `2` for cylindrical, `3` for
  spherical.

- Grid:

  If specified: either an object, as returned by `setup.grid.1D` from
  the package `ReacTran`, a vector of length 101 with the transport
  distances (from mid to mid of layers, upper interface = diffusive
  boundary layer), or one number with the constant layer thickness. If
  `NULL`, it is defined as
  `setup.grid.1D(x.up = 0, dx.1 = 0.01, N = 100, L = 100)`, i.e. the
  total length is 100 cm, the first layer is 0.01 cm thick and layers
  are increasing with depth for 100 layers.

- porosity:

  If specified, either an object with porosities (\[-\]) as returned by
  `setup.prop.1D` from the package `ReacTran`, a vector of length 101
  with the porosities defined at the layer interfaces, or one number
  with the constant porosity. If `NULL`, it is defined by the parameters
  `por0`, `pordeep` and `porcoeff` as:
  `(pordeep + (por0 - pordeep) * exp(-pmax(0, x.int)/porcoeff))`, where
  `x.int` is the distance, from the surface of the layer interface. Note
  that the porosity values should be consistent withe the `Grid` - and
  should be inbetween 0 and 1.

- bioturbation:

  If specified, either an object with bioturbation rates (units
  \[cm2/d\]) as returned by `setup.prop.1D` from the package `ReacTran`,
  a vector of length 101 with the bioturbation defined at the layer
  interfaces, or one number with the constant bioturbation. If `NULL`,
  it is defined by the parameters `biot`, `biotdepth` and `biotatt` as:
  `biot * exp(-pmax(0, (x.int-biotdpeth))/biotatt)`, where `x.int` is
  the distance, from the surface of the layer interface. Note that the
  bioturbation values should be consistent with the `Grid`.

- irrigation:

  If specified, either an object with irrigation rates (units \[/d\]) as
  returned by `setup.prop.1D` from the package `ReacTran`, a vector of
  length 100 with the irrigation rates defined at the layer centres, or
  one number with the constant rates. If `NULL`, it is defined by the
  parameters `irr`, `irrdepth` and `irratt` as:
  `irr * exp(-pmax(0, (x-irrdepth)/irratt))`, where `x` is the distance,
  from the surface, of the layer centres. Note that the irrigation
  values should be consistent with the `Grid`.

- surface:

  If specified, either an object with surface areas (units \[cm2\]) as
  returned by `setup.prop.1D` from the package `ReacTran`, a vector of
  length 101 with the surface areas defined at the layer interfaces, or
  one number with the constant surface area. If `NULL`, it is defined by
  the parameter `gridtype`, and the `Grid` as: `surface = 1` for
  `gridtype == 1`, `surface = rev(2*pi*Grid$x.int)` for `gridtype == 2`
  and `surface = rev(pi*(Grid$x.int)^2)` for `gridtype == 3`.

  Note that the surface values should be consistent with the `Grid`.

- diffusionfactor:

  The multiplication factor necessary to go from molecular diffusion to
  effective sediment diffusion, i.e. that takes into account tortuosity.
  If specified, either an object with these factors (\[-\]) as returned
  by `setup.prop.1D` from the package `ReacTran`, a vector of length 101
  with these factors defined at the layer interfaces, or one number with
  the constant factor. If `NULL`, it is set equal to the porosity. Note
  that the factors should be consistent with the `Grid`.

- yini :

  Initial guess of the steady-state solution.

- dynamicbottomwater:

  If `TRUE`, then the concentrations in the water overlying the sediment
  will also be dynamically described, and with water height equal to
  `Hwater`. Note that this will slow down the simulation.

- ratefactor:

  `NULL` or a list, detailing the forcing function for the
  biogeochemical rate multiplication factor. If not specified (or
  `NULL`), then it is assumed to be 1, if a data series, then the mean
  of the data will be used. If a `list`, it should contain either a data
  time series (`list(data = )`) or parameters determining the
  periodicity of the seasonal signal (defined as
  `list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)`.
  see details.

- method:

  The method to be used for estimating steady-state, will be passed as
  argument to the solver
  [steady.1D](https://rdrr.io/pkg/rootSolve/man/steady.1D.html). The
  default is to use the "stode" method. Other options are to use
  "stodes", "runsteady" or "mix". The last option combines runsteady
  with stode, i.e. the model is first dynamically run for the times
  specified in argument `times`, and then the final value from this run
  used as initial guess for the steady-state estimated using `stode`

- times:

  start and end time of the dynamic run - only if `method` is one of
  `"runsteady"` or `"mix"`.

- verbose:

  If TRUE, will write progession to the screen .

- calcpH:

  if TRUE, the pH will be calculated - note that this will not include
  the effects of calcium carbonate precipitation or dissolution on pH.
  Note also that the pH can be estimated afterwards by running
  [FESDIApH](https://stanleesocca.github.io/FESDIA/reference/FESDIAvars.md).
  See examples.

- ... :

  Any argument passed to the steady-state solver.

## Value

`FESDIAsolve` returns an object of class `FESDIAstd` or `PHDIAstd`, and
of class `steady1D`, as generated by the solvers from R-package
`rootSolve`
([steady.1D](https://rdrr.io/pkg/rootSolve/man/steady.1D.html)\[rootSolve\]).

It contains a.o. the elements:

- `y`, with the state variables at steady-state
  (`FDET, SDET`,`O2, NO3, NH3`, `PO4, FeP, CaP`, `DIC, Fe, FeOH3`,
  `H2S, SO4, CH4, ALK`).

- `O2flux, O2deepflux, NO3flux, NO3deepflux, ...`,
  `FDETflux, FDETdeepflux, SDETflux, SDETdeepflux, FePdeepflux, CaPdeepflux, FeOH3deepflux, OrgCflux, OrgNflux, OrgPflux`,
  the sediment-water and burial fluxes, in nmol/cm2/d.

- `DINDIPflux, DINDIPmean, DINDIPdeep`, the dissolved nitrogen to
  phosphorus ratio of flux, sediment concentrations and deep (burial)
  concentration.

- `TotMin, TotOxic, TotDenit, TotFered TotBSR, TotMeth`, total
  mineralisation, total oxic mineralisation, denitrification, iron
  reduction, biological sulphate reduction, and methanogenesis, in
  nmol/cm2/d.

- `PartOxic, PartDenit, PartFered, PartBSR, PartMethano`, the fraction
  of mineralisation due to oxic, denitrification, iron reduction,
  biological sulphate reduction, and methanogenesis.

- `TotNitri, TotFeoxid, TotH2Soxid, TotCH4oxid, TotAOM, TotFePprod, TotCaPprod, TotFePdesorp, TotCaPdiss, TotNprod, TotPprod, TotNH3ads`,
  integrated rates, nmol/cm2/d.

- `PartPremoved, PartNremoved`, the total P and N removed, relative to
  its production.

- `TOC`, the Total organic carbon concentration profile, %.

- `Cprod, Nprod, Pprod, Oxicmin, Denitrific, Feredmin, BSRmin, Methmin, nitri,Feoxid, H2Soxid, CH4oxid, AOM, FeSprod, FePadsorp, CaPprod`,
  rate profiles, nm/cm3 liquid/d.

- `FePdesorp, CaPdiss`, rate profiles, nm/cm3 solid/d.

The meaning and units of these columns can be assessed via the
R-functions:

[`FESDIAsvar()`](https://stanleesocca.github.io/FESDIA/reference/FESDIAvars.md),
[`FESDIA1D()`](https://stanleesocca.github.io/FESDIA/reference/FESDIAvars.md),
[`FESDIA0D()`](https://stanleesocca.github.io/FESDIA/reference/FESDIAvars.md).
See
[FESDIA0D](https://stanleesocca.github.io/FESDIA/reference/FESDIAvars.md).

## Author

Karline Soetaert

## Details

To solve the model, a steady-state solver from package rootSolve (here
we used [`steady.1D`](https://rdrr.io/pkg/rootSolve/man/steady.1D.html))
is used.

## References

Soetaert K, PMJ Herman and JJ Middelburg, 1996a. A model of early
diagenetic processes from the shelf to abyssal depths. Geochimica
Cosmochimica Acta, 60(6):1019-1040.

Soetaert K, PMJ Herman and JJ Middelburg, 1996b. Dynamic response of
deep-sea sediments to seasonal variation: a model. Limnol. Oceanogr.
41(8): 1651-1668.

## Examples

``` r
#===========================================
# Show parameter values
#===========================================

  FESDIAparms()
#>                      parms         units
#> Cflux         4.566210e+02   nmolC/cm2/d
#> pFast         9.000000e-01             -
#> FeOH3flux     1.000000e+00    nmol/cm2/d
#> CaPflux       0.000000e+00    nmol/cm2/d
#> rFast         6.849315e-02            /d
#> rSlow         1.369863e-04            /d
#> NCrFdet       1.509434e-01     molN/molC
#> NCrSdet       1.509434e-01     molN/molC
#> PCrFdet       9.433962e-03     molP/molC
#> PCrSdet       9.433962e-03     molP/molC
#> BCupLiq       2.000000e+00             -
#> BCdownLiq     3.000000e+00             -
#> O2bw          3.000000e+02       mmol/m3
#> NO3bw         1.000000e+01       mmol/m3
#> NO2bw         0.000000e+00       mmol/m3
#> NH3bw         1.000000e+00       mmol/m3
#> CH4bw         0.000000e+00       mmol/m3
#> PO4bw         5.000000e-01       mmol/m3
#> DICbw         2.100000e+03       mmol/m3
#> Febw          0.000000e+00       mmol/m3
#> H2Sbw         0.000000e+00       mmol/m3
#> SO4bw         3.100000e+04       mmol/m3
#> ALKbw         2.400000e+03       mmol/m3
#> Mnbw          0.000000e+00       mmol/m3
#> O2dw                    NA       mmol/m3
#> NO3dw                   NA       mmol/m3
#> NO2dw                   NA       mmol/m3
#> NH3dw                   NA       mmol/m3
#> CH4dw                   NA       mmol/m3
#> PO4dw                   NA       mmol/m3
#> DICdw                   NA       mmol/m3
#> Fedw                    NA       mmol/m3
#> H2Sdw                   NA       mmol/m3
#> SO4dw                   NA       mmol/m3
#> ALKdw                   NA       mmol/m3
#> Mndw                    NA       mmol/m3
#> w             2.739726e-06          cm/d
#> biot          2.739726e-03         cm2/d
#> biotdepth     5.000000e+00            cm
#> biotatt       1.000000e+00           /cm
#> irr           0.000000e+00            /d
#> irrdepth      5.000000e+00            cm
#> irratt        1.000000e+00            cm
#> gasflux       0.000000e+00          cm/d
#> NH3Ads        1.300000e+00             -
#> rnitri1       2.000000e+01            /d
#> rnitri2       2.000000e+01            /d
#> ranammox      1.000000e-01  /(mmol/m3)/d
#> ksO2nitri     1.000000e+00     mmolO2/m3
#> ksO2oxic      3.000000e+00     mmolO2/m3
#> ksNO3denit    3.000000e+01    mmolNO3/m3
#> kinO2denit    1.000000e+00     mmolO2/m3
#> kinNO3anox    1.000000e+00    mmolNO3/m3
#> kinO2anox     1.000000e-03     mmolO2/m3
#> temperature   1.000000e+01           dgC
#> salinity      3.500000e+01           psu
#> TOC0          5.000000e-01             %
#> rFePadsorp    1.000000e-06            /d
#> rCaPprod      0.000000e+00            /d
#> rCaPdiss      0.000000e+00            /d
#> CPrCaP        2.869565e-01       mol/mol
#> rPads         0.000000e+00            /d
#> rPdes         0.000000e+00            /d
#> maxPads       1.000000e+03 mmolP/m3solid
#> ksFeOH3       1.250000e+04  mmolFeOH3/m3
#> kinFeOH3      1.250000e+04  mmolFeOH3/m3
#> ksSO4BSR      1.600000e+03      mmolS/m3
#> kinSO4Met     1.000000e+03      mmolS/m3
#> rFeox         3.000000e-01  /(mmol/m3)/d
#> rH2Sox        5.000000e-04  /(mmol/m3)/d
#> rFeS          1.000000e-03  /(mmol/m3)/d
#> rCH4ox        2.700000e+01  /(mmol/m3)/d
#> rAOM          3.000000e-05  /(mmol/m3)/d
#> rSurfH2Sox    0.000000e+00            /d
#> rSurfCH4ox    0.000000e+00            /d
#> ksSurfALK     3.000000e+03       mmol/m3
#> ksO2reox      1.000000e+00     mmolO2/m3
#> ODUoxdepth    5.000000e+00            cm
#> ODUoxatt      1.000000e+00           /cm
#> por0          9.000000e-01             -
#> pordeep       5.000000e-01             -
#> porcoeff      3.000000e-01            cm
#> formationtype 1.000000e+00             -
#> dilution      0.000000e+00            /d
#> Hwater        1.000000e+01            cm
#> Cfall         1.000000e+02          cm/d
#> FePfall       1.000000e+02          cm/d
#> FeOH3fall     1.000000e+02          cm/d
#> CaPfall       1.000000e+02          cm/d
#> addalk        1.000000e+00             -
#> MPBprod       0.000000e+00     mmol/m3/d
#> kMPB          4.000000e+00           /cm
#> kDINupt       1.000000e-02       mmol/m3
#> kPO4upt       1.000000e-03       mmol/m3
#> kDICupt       1.000000e+00       mmol/m3
#> rH2Sfeox      1.200000e-04    cm3/nmol/d
#> MnO2flux      1.000000e+00    nmol/cm2/d
#> rAgeFeox      1.555200e-03    cm3/nmol/d
#> rMnOxid       8.640000e-04    cm3/nmol/d
#> rH2SMnox      1.728000e-04    cm3/nmol/d
#> rAgeMnox      4.665600e-03    cm3/nmol/d
#> rMnFe         6.480000e-06    cm3/nmol/d
#> rMnS          1.000000e-05    cm3/nmol/d
#> rMnCO3prec    3.000000e-04    cm3/nmol/d
#> rFeCO3prec    3.000000e-04    cm3/nmol/d
#> ksMnO2        2.600000e+03       mmol/m3
#> pFastFeOx     5.000000e-01             -
#> pFastMnOx     5.000000e-01             -
#> kinMnO2       2.600000e+03       mmol/m3
#> isDICcorr     0.000000e+00             -
#>                                                            description
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

#===========================================
# Runs with different carbon fluxes
#===========================================

  out  <- FESDIAsolve(parms = list(Cflux = 2*1e5/12/365), calcpH = TRUE)
  out2 <- FESDIAsolve(parms = list(Cflux = 20*1e5/12/365), calcpH = TRUE)
  out3 <- FESDIAsolve(parms = list(Cflux = 200*1e5/12/365), calcpH = TRUE)

  par(mar = c(3,3,3,2))
  plot(out, out2, out3, ylim = c(20, 0), mfrow = c(4, 4), 
    which = c(1:14,16), lwd = 2, lty = 1)

  plot(out, out2, out3, which = "pH", ylim = c(20,0))



#===========================================
# Runs that are difficult to solve
#===========================================

# simple run
  out1 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 0))

# simple run used as initial condition for second run
  out2 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 1e3), yini = out1$y)

# second run used as initial condition for third run
# use mixed method: first dynamic run, then steady-state solver
  out3 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 1e4), yini = out2$y, method = "mixed")

  out4 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 5e4), yini = out3$y, 
     method = "runsteady", times = c(0, 1e6))
  
  FESDIAbudgetO2(out1, out2, out3, out4, which = "Rates") 
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#> Warning: number of columns of result, 1, is not a multiple of vector length 2 of arg 2
#>                            [,1]         [,2]         [,3]         [,4]
#> Nitrification      7.321685e+01 8.644882e+01 1.307398e+02 2.027087e+02
#> FeOxidation        3.137167e-01 3.613416e-01 5.655291e-01 7.698778e-01
#> MnOxidation        4.907737e-03 5.225579e-03 6.165785e-03 7.629481e-03
#> H2Soxidation       2.494271e-01 3.190613e-01 1.168650e+00 6.405278e+00
#> CH4oxidation       1.010089e-01 2.623617e-01 3.310563e+00 2.010633e+01
#> H2Soxid.dist       0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> CH4oxid.dist       0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> OxicMineralisation 4.039145e+02 5.159460e+02 1.507571e+03 5.848236e+03
#> MPBO2production    0.000000e+00 1.248485e+02 1.248406e+03 6.240432e+03
#> MPBO2respiration   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> Total              4.778004e+02 7.281913e+02 2.891768e+03 1.231867e+04
```
