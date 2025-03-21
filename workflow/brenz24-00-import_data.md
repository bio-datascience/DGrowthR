brenz24-00-import_data
================
Compiled at 2025-03-21 16:22:06 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "58da8ff6-9fb2-46a3-a33b-89c52f328351")
```

The purpose of this document is …

``` r
library("tidyverse")
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library("DGrowthR")
```

``` r
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

## Read metadata

First we can read and write the metadata data from the brenzinger 2024
study.

``` r
file.copy(from = "../Source_Data/Brenzinger_2024/Database.txt", path_target("Database.txt"))
```

    ## [1] TRUE

``` r
file.copy(from = "../Source_Data/Brenzinger_2024/Map.txt", path_target("Map.txt"))
```

    ## [1] TRUE

## Read the OD data.

``` r
dir.create(path_target("od_data"))

input_dir <- "../Source_Data/Brenzinger_2024/od_data/"
od_files <- list.files(input_dir)

file.copy(file.path(input_dir, od_files), file.path(path_target("od_data"), od_files))
```

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE

## Files written

These files have been written to the target directory,
`data/brenz24-00-import_data`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 3 × 4
    ##   path         type             size modification_time  
    ##   <fs::path>   <fct>     <fs::bytes> <dttm>             
    ## 1 Database.txt file              397 2025-03-21 16:22:10
    ## 2 Map.txt      file            6.88K 2025-03-21 16:22:10
    ## 3 od_data      directory           8 2025-03-21 16:22:10

## Session Info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /opt/bayresq.net/R/R-4.2.0/lib/R/lib/libRblas.so
    ## LAPACK: /opt/bayresq.net/R/R-4.2.0/lib/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] DGrowthR_1.0    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
    ##  [5] dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1    
    ##  [9] tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] mclust_6.1.1        Rcpp_1.0.12         mvtnorm_1.2-4      
    ##  [4] here_1.0.1          fdapace_0.5.9       lattice_0.20-45    
    ##  [7] rprojroot_2.0.4     digest_0.6.35       foreach_1.5.2      
    ## [10] utf8_1.2.4          R6_2.5.1            backports_1.4.1    
    ## [13] evaluate_0.23       pracma_2.4.4        pillar_1.9.0       
    ## [16] rlang_1.1.4         rstudioapi_0.16.0   data.table_1.15.4  
    ## [19] rpart_4.1.16        Matrix_1.6-5        checkmate_2.3.1    
    ## [22] rmarkdown_2.27      foreign_0.8-82      htmlwidgets_1.6.4  
    ## [25] uwot_0.2.2          munsell_0.5.1       compiler_4.2.0     
    ## [28] numDeriv_2016.8-1.1 xfun_0.44           pkgconfig_2.0.3    
    ## [31] base64enc_0.1-3     htmltools_0.5.8.1   nnet_7.3-17        
    ## [34] tidyselect_1.2.1    gridExtra_2.3       htmlTable_2.4.2    
    ## [37] Hmisc_5.1-2         codetools_0.2-18    fansi_1.0.6        
    ## [40] crayon_1.5.2        tzdb_0.4.0          withr_3.0.1        
    ## [43] MASS_7.3-56         grid_4.2.0          jsonlite_1.8.8     
    ## [46] gtable_0.3.5        lifecycle_1.0.4     magrittr_2.0.3     
    ## [49] scales_1.3.0        cli_3.6.3           stringi_1.8.4      
    ## [52] fs_1.6.4            projthis_0.0.0.9025 ellipsis_0.3.2     
    ## [55] generics_0.1.3      vctrs_0.6.5         Formula_1.2-5      
    ## [58] iterators_1.0.14    tools_4.2.0         glue_1.7.0         
    ## [61] hms_1.1.2           fastmap_1.2.0       yaml_2.3.8         
    ## [64] timechange_0.2.0    colorspace_2.1-1    cluster_2.1.3      
    ## [67] knitr_1.46
