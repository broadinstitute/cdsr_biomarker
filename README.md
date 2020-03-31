cds\_biomarker
================

cds\_biomarker is an R toolkit for biomarker analysis. It includes
helpful functions and templates for standard analyses.

## Install

To install cds\_biomarker clone this repository then run this
command:

``` r
install.packages("PATH_TO_LIBRARY/cds_biomarker", repos = NULL, type = "source")
```

## Templates

The
[templates](https://github.com/broadinstitute/cds_biomarker/tree/master/templates)
directory contains templates for standard
    analyses.

  - [multi\_profile\_comparision](https://github.com/broadinstitute/cds_biomarker/tree/master/templates)
    compares biomarkers for multiple response profiles e.g drug and
    genetic or multiple drugs.

## Linear associations

`lin_associations()` can be used to calcualte the linear associations
between a vector of responses and a matrix of biomarker features.

``` r
omics <- load_omics()
A <- omics[["RNA expression"]]
y <- load_achilles("PAX8")
```

``` r
cls <- intersect(names(y),rownames(A))
lin_associations(A[cls,],y[cls]) %>% arrange(p.val) %>% head(5)
```

    ## # A tibble: 5 x 10
    ##   feat.A     p.val  q.val.BH     n betahat sebetahat      lfdr    qvalue
    ##   <chr>      <dbl>     <dbl> <int>   <dbl>     <dbl>     <dbl>     <dbl>
    ## 1 GE_BX… 2.62e-109 1.34e-104   494  -0.230    0.0103 3.64e-105 3.64e-105
    ## 2 GE_CE… 1.51e- 79 3.86e- 75   494  -0.201    0.0106 1.71e- 75 8.57e- 76
    ## 3 GE_PA… 5.53e- 75 9.39e- 71   494  -0.284    0.0155 4.24e- 71 1.41e- 71
    ## 4 GE_KC… 2.71e- 63 3.46e- 59   494  -0.188    0.0112 2.48e- 59 6.19e- 60
    ## 5 GE_RH… 2.38e- 60 2.43e- 56   494  -0.189    0.0115 2.07e- 56 4.14e- 57
    ## # … with 2 more variables: PosteriorMean <dbl>, PosteriorSD <dbl>
