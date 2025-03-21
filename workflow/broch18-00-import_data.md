broch18-00-import_data
================
Compiled at 2025-03-21 17:57:50 UTC

``` r
here::i_am(paste0(params$name, ".Rmd"), uuid = "c1cbddf4-dfaa-4c20-bc79-297c65de4d98")
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
# create or *empty* the target directory, used to write this file's data: 
projthis::proj_create_dir_target(params$name, clean = TRUE)

# function to get path to target directory: path_target("sample.csv")
path_target <- projthis::proj_path_target(params$name)

# function to get path to previous data: path_source("00-import", "sample.csv")
path_source <- projthis::proj_path_source(params$name)
```

``` r
ecoli.data <- read.csv("../Source_Data/Brochado_2018/Ecoli.csv")
```

``` r
ecoli_BW <- ecoli.data %>% filter(Strain == "Ecoli_BW")
ecoli_BW <- select(ecoli_BW, !contains("A22"))

# Pivot longer
ecoli_BW <- ecoli_BW %>%
  # Match query and colnames format
  mutate(Query = gsub(" ", ".", Query)) %>% 
  
  pivot_longer(
    cols = -c(Strain, Query, Time_h),
    names_to = "treatment",
    values_to = "od"
  )

# Prepare OD
ecoli_BW_od <- ecoli_BW %>% 
  unite("curve_id", c(Query, treatment), sep=":") %>% 
  select(-Strain) %>% 
  arrange(curve_id, Time_h) %>% 
  
  group_by(curve_id) %>% 
  mutate(timepoint = seq(from=0, by=0.6, length.out=n()), # Measurements (on average) are done every 
         timepoint_n = 1:n())

# Prepare metadata
ecoli_BW_metadata <- ecoli_BW %>% 
  select(Query, treatment) %>% 
  distinct() %>% 
  
  mutate(query_chem = str_extract(Query, "(.*?)\\.\\d*\\.rep\\.\\d+", group=1),
         query_chem = if_else(is.na(query_chem), str_extract(Query, "(.*?)\\.rep\\.\\d+", group=1), query_chem),
         
         query_conc = str_extract(Query, ".*?\\.(\\d*)\\.rep\\.\\d+", group=1),
         query_rep = str_extract(Query, ".*?\\.\\d*\\.rep\\.(\\d+)", group=1),
         
         treatment_chem = str_extract(treatment, "(.*?)\\.\\d*\\.rep\\.\\d+", group=1),
         treatment_chem = if_else(is.na(treatment_chem), str_extract(treatment, "(.*?)\\.rep\\.\\d+", group=1), treatment_chem),
         
         treatment_conc = str_extract(treatment, ".*?\\.(\\d*)\\.rep\\.\\d+", group=1),
         treatment_rep = str_extract(treatment, ".*?\\.\\d*\\.rep\\.(\\d+)", group=1)) %>% 
  
  rename(query = Query) %>% 
  arrange(query, treatment) %>% 
  
  unite("curve_id", c(query, treatment), sep=":", remove = FALSE)


# Metadata timepoints
tps_by_cid <- ecoli_BW_od %>% 
  group_by(curve_id) %>% 
  summarise(n_tps = n())

# Update metadata
ecoli_BW_metadata <- ecoli_BW_metadata %>% 
  left_join(tps_by_cid, by="curve_id") %>% 
  
  unite("query.combo", c(query_chem, query_conc), sep=".", remove=FALSE) %>% 
  unite("treatment.combo", c(treatment_chem, treatment_conc), sep=".", remove=FALSE) %>% 
  
  unite("combination.id", c(query.combo, treatment.combo), sep = ":", remove=FALSE)
```

``` r
ecoli_iAi1 <- ecoli.data %>% filter(Strain == "Ecoli_iAi1")
ecoli_iAi1 <- select(ecoli_iAi1, !contains("A22"))

# Pivot longer
ecoli_iAi1 <- ecoli_iAi1 %>%
  # Match query and colnames format
  mutate(Query = gsub(" ", ".", Query)) %>% 
  
  pivot_longer(
    cols = -c(Strain, Query, Time_h),
    names_to = "treatment",
    values_to = "od"
  )

# Prepare OD
ecoli_iAi1_od <- ecoli_iAi1 %>% 
  unite("curve_id", c(Query, treatment), sep=":") %>% 
  select(-Strain) %>% 
  arrange(curve_id, Time_h) %>% 
  
  group_by(curve_id) %>% 
  mutate(timepoint = seq(from=0, by=0.6, length.out=n()), # Measurements (on average) are done every 
         timepoint_n = 1:n())

# Prepare metadata
ecoli_iAi1_metadata <- ecoli_iAi1 %>% 
  select(Query, treatment) %>% 
  distinct() %>% 
  
 mutate(query_chem = str_extract(Query, "(.*?)\\.\\d*\\.rep\\.\\d+", group=1),
         query_chem = if_else(is.na(query_chem), str_extract(Query, "(.*?)\\.rep\\.\\d+", group=1), query_chem),
         
         query_conc = str_extract(Query, ".*?\\.(\\d*)\\.rep\\.\\d+", group=1),
         query_rep = str_extract(Query, ".*?\\.\\d*\\.rep\\.(\\d+)", group=1),
         
         treatment_chem = str_extract(treatment, "(.*?)\\.\\d*\\.rep\\.\\d+", group=1),
         treatment_chem = if_else(is.na(treatment_chem), str_extract(treatment, "(.*?)\\.rep\\.\\d+", group=1), treatment_chem),
         
         treatment_conc = str_extract(treatment, ".*?\\.(\\d*)\\.rep\\.\\d+", group=1),
         treatment_rep = str_extract(treatment, ".*?\\.\\d*\\.rep\\.(\\d+)", group=1)) %>% 
  
  rename(query = Query) %>% 
  arrange(query, treatment) %>% 
  
  unite("curve_id", c(query, treatment), sep=":", remove = FALSE)


# Metadata timepoints
tps_by_cid <- ecoli_iAi1_od %>% 
  group_by(curve_id) %>% 
  summarise(n_tps = n())

# Update metadata
ecoli_iAi1_metadata <- ecoli_iAi1_metadata %>% 
  left_join(tps_by_cid, by="curve_id") %>% 
  
  unite("query.combo", c(query_chem, query_conc), sep=".", remove=FALSE) %>% 
  unite("treatment.combo", c(treatment_chem, treatment_conc), sep=".", remove=FALSE) %>% 
  
  unite("combination.id", c(query.combo, treatment.combo), sep = ":", remove=FALSE)
```

## Write files

``` r
write_tsv(ecoli_BW_od, path_target("ecoli_od_BW.tsv.gz"))
write_tsv(ecoli_BW_metadata, path_target("ecoli_metadata_BW.tsv.gz"))

write_tsv(ecoli_iAi1_od, path_target("ecoli_od_iAi.tsv.gz"))
write_tsv(ecoli_iAi1_metadata, path_target("ecoli_metadata_iAi.tsv.gz"))
```

## Files written

These files have been written to the target directory,
`data/broch18-00-import_data`:

``` r
projthis::proj_dir_info(path_target())
```

    ## # A tibble: 4 × 4
    ##   path                      type         size modification_time  
    ##   <fs::path>                <fct> <fs::bytes> <dttm>             
    ## 1 ecoli_metadata_BW.tsv.gz  file        1.33M 2025-03-21 17:58:25
    ## 2 ecoli_metadata_iAi.tsv.gz file        1.42M 2025-03-21 17:58:41
    ## 3 ecoli_od_BW.tsv.gz        file       19.37M 2025-03-21 17:58:19
    ## 4 ecoli_od_iAi.tsv.gz       file       20.78M 2025-03-21 17:58:34

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
    ##  [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
    ##  [5] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
    ##  [9] ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] pillar_1.9.0        compiler_4.2.0      tools_4.2.0        
    ##  [4] bit_4.0.5           digest_0.6.35       timechange_0.2.0   
    ##  [7] evaluate_0.23       lifecycle_1.0.4     gtable_0.3.5       
    ## [10] pkgconfig_2.0.3     rlang_1.1.4         cli_3.6.3          
    ## [13] rstudioapi_0.16.0   parallel_4.2.0      yaml_2.3.8         
    ## [16] xfun_0.44           fastmap_1.2.0       withr_3.0.1        
    ## [19] knitr_1.46          fs_1.6.4            projthis_0.0.0.9025
    ## [22] generics_0.1.3      vctrs_0.6.5         hms_1.1.2          
    ## [25] bit64_4.0.5         rprojroot_2.0.4     grid_4.2.0         
    ## [28] tidyselect_1.2.1    glue_1.7.0          here_1.0.1         
    ## [31] R6_2.5.1            fansi_1.0.6         vroom_1.6.1        
    ## [34] rmarkdown_2.27      tzdb_0.4.0          magrittr_2.0.3     
    ## [37] scales_1.3.0        htmltools_0.5.8.1   ellipsis_0.3.2     
    ## [40] colorspace_2.1-1    utf8_1.2.4          stringi_1.8.4      
    ## [43] munsell_0.5.1       crayon_1.5.2
