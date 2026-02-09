# Model solutions to aid in solving new applications

`FESDIAdefault, FESDIAdefaults2` are two different FESDIA model
solutions.

## Usage

``` r
data(FESDIAdefault)
```

## References

Soetaert K, PMJ Herman and JJ Middelburg, 1996a. A model of early
diagenetic processes from the shelf to abyssal depths. Geochimica
Cosmochimica Acta, 60(6):1019-1040.

Soetaert K, PMJ Herman and JJ Middelburg, 1996b. Dynamic response of
deep-sea sediments to seasonal variation: a model. Limnol. Oceanogr.
41(8): 1651-1668.

## Examples

``` r
# defaults
  head(FESDIAparms(FESDIAdefault, as.vector = TRUE)) 
#>        Cflux        pFast    FeOH3flux      CaPflux        rFast        rSlow 
#> 4.566210e+02 9.000000e-01 1.000000e+00 0.000000e+00 6.849315e-02 1.369863e-04 
#        FESDIAparms(FESDIAdefault2))
#  plot(FESDIAdefault, xyswap = TRUE, ylim = c(30, 0),
#    mfrow = c(5,4), grid = FESDIAdepth(FESDIAdefault))
  
```
