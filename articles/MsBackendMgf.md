# Description and usage of MsBackendMgf

**Package**:
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)*\
**Authors**: RforMassSpectrometry Package Maintainer \[cre\], Laurent
Gatto \[aut\] (ORCID: <https://orcid.org/0000-0002-1520-2268>), Johannes
Rainer \[aut\] (ORCID: <https://orcid.org/0000-0002-6977-7147>),
Sebastian Gibb \[aut\] (ORCID: <https://orcid.org/0000-0001-7406-4443>),
Michael Witting \[ctb\] (ORCID:
<https://orcid.org/0000-0002-1462-4426>), Adriano Rutz \[ctb\] (ORCID:
<https://orcid.org/0000-0003-0443-9902>), Corey Broeckling \[ctb\]
(ORCID: <https://orcid.org/0000-0002-6158-827X>)\
**Last modified:** 2026-03-02 14:41:42.922087\
**Compiled**: Mon Mar 2 14:46:26 2026

## Introduction

The *[Spectra](https://bioconductor.org/packages/3.23/Spectra)* package
provides a central infrastructure for the handling of Mass Spectrometry
(MS) data. The package supports interchangeable use of different
*backends* to import MS data from a variety of sources (such as mzML
files). The
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)*
package allows the import of MS/MS data from MGF ([Mascot Generic
Format](http://www.matrixscience.com/help/data_file_help.md)) files. The
`MsBackendMgf` backend allows to load and represent data from these
*standard* MGF files, while the `MsBackendAnnotatedMgf` backend supports
import of data from MGF files containing also annotations for the
individual mass peaks. This vignette illustrates the usage of the
*MsBackendMgf* package.

## Installation

To install this package, start `R` and enter:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MsBackendMgf")
```

This will install this package and all eventually missing dependencies.

## Importing MS/MS data from MGF files

MGF files store one to multiple spectra, typically centroided and of MS
level 2. In our short example below, we load MGF files which are
provided with this package. In a next section we also import data from
an MGF file containing individual peak annotations. Below we first load
all required packages and define the paths to the MGF files.

``` r

library(Spectra)
library(MsBackendMgf)

fls <- dir(system.file("extdata", package = "MsBackendMgf"),
           full.names = TRUE, pattern = "^spectra(.*).mgf$")
fls
```

    ## [1] "/__w/_temp/Library/MsBackendMgf/extdata/spectra.mgf"             
    ## [2] "/__w/_temp/Library/MsBackendMgf/extdata/spectra2.mgf"            
    ## [3] "/__w/_temp/Library/MsBackendMgf/extdata/spectra3_empty_peaks.mgf"
    ## [4] "/__w/_temp/Library/MsBackendMgf/extdata/spectra4.mgf"

MS data can be accessed and analyzed through `Spectra` objects. Below we
create a `Spectra` with the data from these MGF files. To this end we
provide the file names and specify to use a
[`MsBackendMgf()`](https://rformassspectrometry.github.io/MsBackendMgf/reference/MsBackendMgf.md)
backend as *source* to enable data import.

``` r

sps <- Spectra(fls, source = MsBackendMgf())
```

    ## Start data import from 4 files ... done

With that we have now full access to all imported spectra variables that
we list below.

``` r

spectraVariables(sps)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "TITLE"                  
    ## [19] "RAWFILE"                 "CLUSTER_ID"             
    ## [21] "MSLEVEL"

Besides default spectra variables, such as `msLevel`, `rtime`,
`precursorMz`, we also have additional spectra variables such as the
`TITLE` of each spectrum in the MGF file.

``` r

sps$rtime
```

    ##  [1] 1028.000 1117.000 1127.000 2678.940 2373.511 2511.030       NA  162.070
    ##  [9] 1028.000 1028.000

``` r

sps$TITLE
```

    ##  [1] "File193 Spectrum1719 scans: 2162"                                                         
    ##  [2] "File193 Spectrum1944 scans: 2406"                                                         
    ##  [3] "File193 Spectrum1968 scans: 2432"                                                         
    ##  [4] "mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840:WNQLQAFWGTGK/2"
    ##  [5] "mzspec:PXD002084:TCGA-AA-A01D-01A-23_W_VU_20121106_A0218_5I_R_FR15:scan:5209:DLTDYLMK/2"  
    ##  [6] "mzspec:MSV000080679:j11962_C1orf144:scan:10671:DLTDYLMK/2"                                
    ##  [7] "CCMSLIB00000840351"                                                                       
    ##  [8] "blank_2-A,1_01_29559.812.812.1"                                                           
    ##  [9] "File193 Spectrum1719 scans: 2162"                                                         
    ## [10] "File193 Spectrum1719 scans: 2162"

By default, fields in the MGF file are mapped to spectra variable names
using the mapping returned by the
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
function:

``` r

spectraVariableMapping(MsBackendMgf())
```

    ##              rtime     acquisitionNum        precursorMz precursorIntensity 
    ##      "RTINSECONDS"            "SCANS"          "PEPMASS"       "PEPMASSINT" 
    ##    precursorCharge 
    ##           "CHARGE"

The names of this `character` vector are the spectra variable names
(such as `"rtime"`) and the field in the MGF file that contains that
information are the values (such as `"RTINSECONDS"`). Note that it is
also possible to overwrite this mapping (e.g. for certain MGF
*dialects*) or to add additional mappings. Below we add the mapping of
the MGF field `"TITLE"` to a spectra variable called `"spectrumName"`.

``` r

map <- c(spectrumName = "TITLE", spectraVariableMapping(MsBackendMgf()))
map
```

    ##       spectrumName              rtime     acquisitionNum        precursorMz 
    ##            "TITLE"      "RTINSECONDS"            "SCANS"          "PEPMASS" 
    ## precursorIntensity    precursorCharge 
    ##       "PEPMASSINT"           "CHARGE"

We can then pass this mapping to the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method, or the
[`Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) constructor.

``` r

sps <- Spectra(fls, source = MsBackendMgf(), mapping = map)
```

    ## Start data import from 4 files ... done

We can now access the spectrum’s title with the newly created spectra
variable `"spectrumName"`:

``` r

sps$spectrumName
```

    ##  [1] "File193 Spectrum1719 scans: 2162"                                                         
    ##  [2] "File193 Spectrum1944 scans: 2406"                                                         
    ##  [3] "File193 Spectrum1968 scans: 2432"                                                         
    ##  [4] "mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840:WNQLQAFWGTGK/2"
    ##  [5] "mzspec:PXD002084:TCGA-AA-A01D-01A-23_W_VU_20121106_A0218_5I_R_FR15:scan:5209:DLTDYLMK/2"  
    ##  [6] "mzspec:MSV000080679:j11962_C1orf144:scan:10671:DLTDYLMK/2"                                
    ##  [7] "CCMSLIB00000840351"                                                                       
    ##  [8] "blank_2-A,1_01_29559.812.812.1"                                                           
    ##  [9] "File193 Spectrum1719 scans: 2162"                                                         
    ## [10] "File193 Spectrum1719 scans: 2162"

In addition we can also access the m/z and intensity values of each
spectrum.

``` r

mz(sps)
```

    ## NumericList of length 10
    ## [[1]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666
    ## [[2]] 101.07074 102.05486 103.00227 ... 1331.56726 1348.58496 1349.59241
    ## [[3]] 102.05556 103.00014 115.05058 ... 1333.599 1334.61304 1335.64368
    ## [[4]] 101.07122 109.68925 115.86999 120.0811 ... 1260.6073 1261.614 1272.6572
    ## [[5]] 130.164459228516 144.150299072266 ... 1019.23852539062 1020.52404785156
    ## [[6]] 110.070594787598 120.080627441406 ... 887.756652832031 998.447387695312
    ## [[7]] 51.022404 57.033543 57.060638 ... 636.130188 660.481445 753.358521
    ## [[8]] numeric(0)
    ## [[9]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666
    ## [[10]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666

``` r

intensity(sps)
```

    ## NumericList of length 10
    ## [[1]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175
    ## [[2]] 1228.93 1424.66 1550.9 1455.45 ... 7380.41 4960.92 5743.83 1780.76
    ## [[3]] 1340.44 1714.76 1938.82 1450.36 2019 ... 5323.02 2265.43 4768.14 1532.12
    ## [[4]] 81011.57 4123.349 4006.9321 66933.17 ... 22042.248 18096.48 12666.438
    ## [[5]] 14.1766004562378 18.5806427001953 ... 22.7096385955811 14.864013671875
    ## [[6]] 1748.57495117188 8689.9951171875 ... 2907.08422851562 2663.30908203125
    ## [[7]] 65.219513 178.758606 13.01786 119.898499 ... 22.05921 30.57095 14.11111
    ## [[8]] numeric(0)
    ## [[9]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175
    ## [[10]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175

The `MsBackendMgf` backend allows also to export data in MGF format.
Below we export the data to a temporary file. We hence call the
[`export()`](https://rdrr.io/pkg/Spectra/man/Spectra.html) function on
our `Spectra` object specifying `backend = MsBackendMgf()` to use this
backend for the export of the data. Note that we use again our custom
mapping of variables such that the spectra variable `"spectrumName"`
will be exported as the spectrums’ title.

``` r

fl <- tempfile()
export(sps, backend = MsBackendMgf(), file = fl, mapping = map)
```

We next read the first lines from the exported file to verify that the
title was exported properly.

``` r

readLines(fl)[1:12]
```

    ##  [1] "BEGIN IONS"                            
    ##  [2] "msLevel=2"                             
    ##  [3] "RTINSECONDS=1028"                      
    ##  [4] "SCANS=2162"                            
    ##  [5] "centroided=TRUE"                       
    ##  [6] "PEPMASS=816.33826"                     
    ##  [7] "CHARGE=2+"                             
    ##  [8] "TITLE=File193 Spectrum1719 scans: 2162"
    ##  [9] "102.0548 753.738"                      
    ## [10] "103.00494 385.376"                     
    ## [11] "103.03531 315.441"                     
    ## [12] "115.05001 413.206"

Note that the `MsBackendMgf` exports all spectra variables as fields in
the mgf file. To illustrate this we add below a new spectra variable to
the object and export the data.

``` r

sps$new_variable <- "A"
export(sps, backend = MsBackendMgf(), file = fl)
readLines(fl)[1:12]
```

    ##  [1] "BEGIN IONS"                                   
    ##  [2] "TITLE=msLevel 2; retentionTime ; scanNum "    
    ##  [3] "msLevel=2"                                    
    ##  [4] "RTINSECONDS=1028"                             
    ##  [5] "SCANS=2162"                                   
    ##  [6] "centroided=TRUE"                              
    ##  [7] "PEPMASS=816.33826"                            
    ##  [8] "CHARGE=2+"                                    
    ##  [9] "spectrumName=File193 Spectrum1719 scans: 2162"
    ## [10] "new_variable=A"                               
    ## [11] "102.0548 753.738"                             
    ## [12] "103.00494 385.376"

We can see that also our newly defined variable was exported. Also,
because we did not provide our custom variable mapping this time, the
variable `"spectrumName"` was **not** used as the spectrum’s title.

Sometimes it might be required to not export all spectra variables since
some exported fields might not be recognized/supported by external
tools. Using the
[`selectSpectraVariables()`](https://rdrr.io/pkg/Spectra/man/filterMsLevel.html)
function we can reduce our `Spectra` object to export to contain only
relevant spectra variables. Below we restrict the data to only m/z,
intensity, retention time, acquisition number, precursor m/z and
precursor charge and export these to an MGF file. Also, some external
tools don’t support the `"TITLE"` field in the MGF file. To disable
export of the spectrum ID/title `exportTitle = FALSE` can be used.

``` r

sps_ex <- selectSpectraVariables(sps, c("mz", "intensity", "rtime",
                                        "acquisitionNum", "precursorMz",
                                        "precursorCharge"))
export(sps_ex, backend = MsBackendMgf(), file = fl, exportTitle = FALSE)
readLines(fl)[1:12]
```

    ##  [1] "BEGIN IONS"        "RTINSECONDS=1028"  "SCANS=2162"       
    ##  [4] "PEPMASS=816.33826" "CHARGE=2+"         "102.0548 753.738" 
    ##  [7] "103.00494 385.376" "103.03531 315.441" "115.05001 413.206"
    ## [10] "115.08686 588.273" "120.08063 800.016" "124.10555 526.761"

### Annotated MGF files

Variants of MGF files can also contain annotations for the individual
mass peaks. These files are expected to contain one line per mass peak,
with the first two elements being the mass peak’s *m/z* and intensity
values. All (eventually present) additional elements are considered
*annotations*. The `MsBackendAnnotatedMgf` backend can be used to import
and represent such files. Below we read the first 10 lines of an example
annotated MGF file.

``` r

fl <- dir(system.file("extdata", package = "MsBackendMgf"),
          pattern = "fiora", full.names = TRUE)
readLines(fl, 10)
```

    ##  [1] "BEGIN IONS"                                                                                                                
    ##  [2] "TITLE=PC(PGF1alpha/22:6(4Z,7Z,10Z,13Z,16Z,19Z))"                                                                           
    ##  [3] "SMILES=CCCCCC@@HO"                                                                                                         
    ##  [4] "PRECURSORTYPE=[M+H]+"                                                                                                      
    ##  [5] "COLLISIONENERGY=10.0"                                                                                                      
    ##  [6] "INSTRUMENTTYPE=HCD"                                                                                                        
    ##  [7] "COMMENT=\"In silico generated spectrum by Fiora OS v0.1.0\""                                                               
    ##  [8] "906.5854755900899 0.24612045288085938 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"
    ##  [9] "57.06987670808999 0.0012485487386584282 CCCC//[M-H]+"                                                                      
    ## [10] "71.08552677208999 0.0015523125184699893 CCCCC//[M-H]+"

The last 3 displayed lines show the mass peak content of a MGF spectrum
with the *m/z*, intensity and an annotation for the mass peak. Such data
files can be imported using the `MsBackendAnntotatedMgf` backend:

``` r

sps_ann <- Spectra(fl, source = MsBackendAnnotatedMgf())
```

    ## Start data import from 1 files ... done

``` r

sps_ann
```

    ## MSn data (Spectra) with 1 spectra in a MsBackendAnnotatedMgf backend:
    ##     msLevel     rtime scanIndex
    ##   <integer> <numeric> <integer>
    ## 1         2        NA        NA
    ##  ... 23 more variables/columns.

The peaks annotation have also been imported and are available as
additional *peaks variables*:

``` r

peaksVariables(sps_ann)
```

    ## [1] "mz"        "intensity" "V1"

The names of the annotations follow the standard R convention of naming
columns in a `data.frame` if no name is provided, i.e. it consists of
`"V"` followed by a number. The data can be retrieved with the
[`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html)
function specifying the names of the peaks variables to extract. If not
provided, `peakData()` will only extract the *m/z* and intensity values.
Thus, below we call
[`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html)
providing the names of all requested variables (columns).

``` r

pd <- peaksData(sps_ann, columns = c("mz", "intensity", "V1"))
```

The peaks data of one spectrum is now represented as a `data.frame`:

``` r

pd[[1L]] |>
    head()
```

    ##         mz   intensity
    ## 1 15.99546 0.004359378
    ## 2 19.01784 0.009557842
    ## 3 53.03858 0.001721356
    ## 4 55.05423 0.002519508
    ## 5 57.06988 0.001248549
    ## 6 57.06988 0.002908136
    ##                                                                                     V1
    ## 1 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+
    ## 2                                                                         CCCC//[M-H]+
    ## 3                                                                        CCCCC//[M-H]+
    ## 4                                                                   CCCCC[CH]O//[M-H]+
    ## 5                                                   C=C[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+
    ## 6   CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M]+

## Parallel processing

The *MsBackendMgf* package supports parallel processing for data import.
Parallel processing can be enabled by providing the parallel processing
setup to the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function (or the
[`readMgf()`](https://rformassspectrometry.github.io/MsBackendMgf/reference/readMgf.md)
function) with the `BPPARAM` parameter. By default (with
`BPPARAM = SerialParam()`) parallel processing is disabled. If enabled,
and data import is performed on a single file, the extraction of spectra
information on the imported MGF file is performed in parallel. If data
is to be imported from multiple files, the import is performed in
parallel on a per-file basis (i.e. in parallel from the different
files). Generally, the performance gain through parallel processing is
only moderate and it is only suggested if a large number of files need
to be processed, or if the MGF file is very large (e.g. containing over
100,000 spectra).

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2026-03-01 r89508)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] MsBackendMgf_1.19.0 Spectra_1.21.2      BiocParallel_1.45.0
    ## [4] S4Vectors_0.49.0    BiocGenerics_0.57.0 generics_0.1.4     
    ## [7] BiocStyle_2.39.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0         compiler_4.6.0         BiocManager_1.30.27   
    ##  [4] parallel_4.6.0         cluster_2.1.8.2        jquerylib_0.1.4       
    ##  [7] systemfonts_1.3.1      IRanges_2.45.0         textshaping_1.0.4     
    ## [10] yaml_2.3.12            fastmap_1.2.0          R6_2.6.1              
    ## [13] ProtGenerics_1.39.2    knitr_1.51             htmlwidgets_1.6.4     
    ## [16] MASS_7.3-65            bookdown_0.46          desc_1.4.3            
    ## [19] bslib_0.10.0           rlang_1.1.7            cachem_1.1.0          
    ## [22] xfun_0.56              fs_1.6.6               MsCoreUtils_1.23.2    
    ## [25] sass_0.4.10            otel_0.2.0             cli_3.6.5             
    ## [28] pkgdown_2.2.0.9000     digest_0.6.39          MetaboCoreUtils_1.19.2
    ## [31] lifecycle_1.0.5        clue_0.3-67            evaluate_1.0.5        
    ## [34] codetools_0.2-20       ragg_1.5.0             rmarkdown_2.30        
    ## [37] tools_4.6.0            htmltools_0.5.9
