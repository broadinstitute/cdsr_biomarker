cdsrbiomarker
================

cdsrbiomarker is an R toolkit for biomarker analysis. It includes
helpful functions and standard reports.

## Install

``` r
library(devtools)
devtools::install_github("broadinstitute/cdsr_models")
devtools::install_github("broadinstitute/cdsr_biomarker")
```

The package can then be loaded by calling

``` r
library(cdsrbiomarker)
```

## Taiga

Many of the functions in the package use data which is stored on taiga.
If you are at the Broad you can install taigr the taiga client for R by
following the instructions
[here](https://github.com/broadinstitute/taigr).

## Functions

### get\_biomarkers

`get_biomarkers` is a wrapper around the `discrete_test`,
`lin_associations`, and `random_forest` functions in `cdsrmodels`. It
takes a response vector `y` and then loads features matricies from
Taiga. The appropriate modeling functions are applied to each feature
matrix and the results are returned and saved in the `out_path`
directory if it is
specified.

``` r
biomarker_results <- get_biomarkers(y, p_cutoff=0.1, out_path="~/Desktop/example")
```

## Reports

The
[reports](https://github.com/broadinstitute/cdsr_biomarker/tree/master/inst/reports)
directory contains standard biomarker
    reports.

  - [single\_profile\_biomarker\_report](https://github.com/broadinstitute/cdsr_biomarker/tree/master/inst/reports/single_profile_biomarker_report.Rmd)
    generates a biomarker report for a single response
    profile.

  - [multi\_profile\_biomarker\_report](https://github.com/broadinstitute/cdsr_biomarker/tree/master/inst/reports/multi_profile_biomarker_report.Rmd)
    compares biomarkers for multiple response profiles e.g drug and
    genetic or multiple drugs.

There are wrapper functions in cdsrbiomarker to automaticaly genenarate
these reports. Here is an example using Achilles data for EGFR and PRISM
data for a few EGFR inhibitors:

1.  Make a cell line by perturbation response matrix
Y.

<!-- end list -->

``` r
gene_effect <- load.from.taiga(data.name='depmap-a0ab',data.file='Achilles_gene_effect')[,"EGFR (1956)"] %>% 
  enframe(name = "arxspan_id",value = "xpr_egfr")
auc <- load.from.taiga(data.name='secondary-screen-15e6', data.file='secondary_merged_drc_parameters') %>% 
  filter(repurposing_name %in% c("erlotinib","gefitinib","lapatinib")) %>% 
  select(auc,repurposing_name,arxspan_id) %>% 
  spread(key = "repurposing_name",value = "auc")
Y <- full_join(gene_effect,auc, by = "arxspan_id") %>% column_to_rownames(var = "arxspan_id") %>% as.matrix()
```

``` r
corner(Y)
```

    ##               xpr_egfr erlotinib gefitinib lapatinib
    ## ACH-000004  0.16303435        NA        NA        NA
    ## ACH-000005  0.05234623        NA        NA        NA
    ## ACH-000007 -0.27242329 1.2249556 0.6756661 0.6937505
    ## ACH-000009 -0.78660538        NA        NA        NA
    ## ACH-000011 -0.10197377 0.7143025 0.7616728 0.7627777

2.  Make a meta data table which will be displayed in the
report.

<!-- end list -->

``` r
meta_data <- list(perturbation = colnames(Y), type = c("CRISPR","Drug","Drug","Drug")) %>% as_tibble()
meta_data
```

    ## # A tibble: 4 x 2
    ##   perturbation type  
    ##   <chr>        <chr> 
    ## 1 xpr_egfr     CRISPR
    ## 2 erlotinib    Drug  
    ## 3 gefitinib    Drug  
    ## 4 lapatinib    Drug

3.  Call the generate report function for the report you want. You will
    need to give it a file path to save the results
to.

<!-- end list -->

``` r
cdsrbiomarker::generate_multi_profile_biomarker_report("~/Desktop/example/","example_title",Y,meta_data)
```

4.  If you already have the biomarker results files and just want to
    generate the report you can do it like
this.

<!-- end list -->

``` r
cdsrbiomarker::generate_multi_profile_biomarker_report("~/Desktop/example","example")
```
