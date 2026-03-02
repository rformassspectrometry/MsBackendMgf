# MS data backend for MGF files

The `MsBackendMgf` class supports import and export of MS/MS spectra
data from/to files in Mascot Generic Format
([mgf](http://www.matrixscience.com/help/data_file_help.md)) files.
After initial import, the full MS data is kept in memory. `MsBackendMgf`
extends the
[`Spectra::MsBackendDataFrame()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
backend directly and supports thus the
[`Spectra::applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
function to make data manipulations persistent.

The `MsBackendAnnotatedMgf` class supports import of data from MGF files
that provide, in addition to the *m/z* and intensity values, also
additional annotations/metadata for each mass peak. For such MGF files
it is expected that each line contains information from a single mass
peak, separated by a white space (blank). The first two elements are
expected to be the peak's *m/z* and intensity values, while each
additional element is considered an annotation for this specific peak.
See examples below for the format of a supported MGF file. The
`backendInitialize()` method of `MsBackendAnnotatedMgf` does not support
parameter `nlines`. Also, import of data can be considerably slower
compared to the standard `MsBackendMgf` backend, because of the
additionally required parsing of peak annotations. Peaks information in
MGF files are not named, thus, additional peaks annotations are named
using the standard naming convention for column named of data frames:
the first peaks annotation is called `"V1"`, the second (if available)
`"V2"` and so on.

New objects are created with the `MsBackendMgf()` or
`MsBackendAnnotatedMgf()` function. The `backendInitialize()` method has
to be subsequently called to initialize the object and import the MS/MS
data from (one or more) MGF files.

The `MsBackendMgf` backend provides an `export` method that allows to
export the data from the `Spectra` object (parameter `x`) to a file in
mgf format. See the package vignette for details and examples.

Default mappings from fields in the MGF file to spectra variable names
are provided by the `spectraVariableMapping()` function. This function
returns a named character vector were names are the spectra variable
names and the values the respective field names in the MGF files. This
named character vector is submitted to the import and export function
with parameter `mapping`. It is also possible to pass own mappings (e.g.
for special MGF dialects) with the `mapping` parameter.

## Usage

``` r
# S4 method for class 'MsBackendMgf'
backendInitialize(
  object,
  files,
  mapping = spectraVariableMapping(object),
  nlines = -1L,
  ...,
  BPPARAM = SerialParam()
)

MsBackendMgf()

# S4 method for class 'MsBackendMgf'
spectraVariableMapping(object, format = c("mgf"))

# S4 method for class 'MsBackendMgf'
export(
  object,
  x,
  file = tempfile(),
  mapping = spectraVariableMapping(object),
  exportTitle = TRUE,
  ...
)

# S4 method for class 'MsBackendAnnotatedMgf'
backendInitialize(
  object,
  files,
  mapping = spectraVariableMapping(object),
  ...,
  BPPARAM = SerialParam()
)

MsBackendAnnotatedMgf()
```

## Arguments

- object:

  Instance of `MsBackendMgf` class.

- files:

  `character` with the (full) file name(s) of the mgf file(s) from which
  MS/MS data should be imported.

- mapping:

  for `backendInitialize()` and `export`: named `character` vector
  allowing to specify how fields from the MGF file should be renamed.
  Names are supposed to be the spectra variable name and values of the
  vector the field names in the MGF file. See output of
  `spectraVariableMapping()` for the expected format and examples below
  or description above for details.

- nlines:

  for `backendInitialize()` of `MsBackendMgf`: `integer(1)` defining the
  number of lines that should be imported and processed from the MGF
  file(s). By default (`nlines = -1L`) the full file is imported and
  processed at once. If set to a positive integer, the data is imported
  and processed *chunk-wise* using
  [`readMgfSplit()`](https://rformassspectrometry.github.io/MsBackendMgf/reference/readMgf.md).

- ...:

  Currently ignored.

- BPPARAM:

  Parameter object defining the parallel processing setup. If parallel
  processing is enabled (with `BPPARAM` different than
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html),
  the default) and length of `files` is larger than one, import is
  performed in parallel on a per-file basis. If data is to be imported
  from a single file (i.e., length of `files` is one), parsing of the
  imported file is performed in parallel. See also
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  for information on available parallel processing setup options.

- format:

  for `spectraVariableMapping()`: `character(1)` defining the format to
  be used. Currently only `format = "mgf"` is supported.

- x:

  for `export()`: an instance of
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  class with the data that should be exported.

- file:

  `character(1)` with the (full) file name to which the data should be
  exported.

- exportTitle:

  `logical(1)` whether the *TITLE* field should be included in the
  exported MGF file. If `TRUE` (the default) a `spectraVariable` called
  `"TITLE"` will be used, if no such variable is present either the
  `spectraNames(object)` will be used or, if they are empty, a title
  will be generated including the MS level, retention time and
  acquisition number of the spectrum.

## Value

See description above.

## Author

Laurent Gatto, Corey Broeckling and Johannes Rainer

## Examples

``` r

library(BiocParallel)
#' Getting the file names of all example MGF files from MsBackendMgf
fls <- dir(system.file("extdata", package = "MsBackendMgf"),
    full.names = TRUE, pattern = "^spectra(.*).mgf$")

## Create an MsBackendMgf backend and import data from test mgf files.
be <- backendInitialize(MsBackendMgf(), fls)
#> Start data import from 4 files ... 
#> done
be
#> MsBackendMgf with 10 spectra
#>      msLevel     rtime scanIndex
#>    <integer> <numeric> <integer>
#> 1          2   1028.00        NA
#> 2          2   1117.00        NA
#> 3          2   1127.00        NA
#> 4          2   2678.94        NA
#> 5          2   2373.51        NA
#> 6          2   2511.03        NA
#> 7          2        NA        NA
#> 8          2    162.07        NA
#> 9          2   1028.00        NA
#> 10         2   1028.00        NA
#>  ... 20 more variables/columns.

be$msLevel
#>  [1] 2 2 2 2 2 2 2 2 2 2
be$intensity
#> NumericList of length 10
#> [[1]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175
#> [[2]] 1228.93 1424.66 1550.9 1455.45 ... 7380.41 4960.92 5743.83 1780.76
#> [[3]] 1340.44 1714.76 1938.82 1450.36 2019 ... 5323.02 2265.43 4768.14 1532.12
#> [[4]] 81011.57 4123.349 4006.9321 66933.17 ... 22042.248 18096.48 12666.438
#> [[5]] 14.1766004562378 18.5806427001953 ... 22.7096385955811 14.864013671875
#> [[6]] 1748.57495117188 8689.9951171875 ... 2907.08422851562 2663.30908203125
#> [[7]] 65.219513 178.758606 13.01786 119.898499 ... 22.05921 30.57095 14.11111
#> [[8]] numeric(0)
#> [[9]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175
#> [[10]] 753.738 385.376 315.441 413.206 ... 3038.73 2016.43 1146.04 704.175
be$mz
#> NumericList of length 10
#> [[1]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666
#> [[2]] 101.07074 102.05486 103.00227 ... 1331.56726 1348.58496 1349.59241
#> [[3]] 102.05556 103.00014 115.05058 ... 1333.599 1334.61304 1335.64368
#> [[4]] 101.07122 109.68925 115.86999 120.0811 ... 1260.6073 1261.614 1272.6572
#> [[5]] 130.164459228516 144.150299072266 ... 1019.23852539062 1020.52404785156
#> [[6]] 110.070594787598 120.080627441406 ... 887.756652832031 998.447387695312
#> [[7]] 51.022404 57.033543 57.060638 ... 636.130188 660.481445 753.358521
#> [[8]] numeric(0)
#> [[9]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666
#> [[10]] 102.0548 103.00494 103.03531 ... 1388.58691 1405.59729 1406.57666

## The spectra variables that are available; note that not all of them
## have been imported from the MGF files.
spectraVariables(be)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "mz"                      "intensity"              
#>  [7] "dataStorage"             "dataOrigin"             
#>  [9] "centroided"              "smoothed"               
#> [11] "polarity"                "precScanNum"            
#> [13] "precursorMz"             "precursorIntensity"     
#> [15] "precursorCharge"         "collisionEnergy"        
#> [17] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [19] "isolationWindowUpperMz"  "TITLE"                  
#> [21] "RAWFILE"                 "CLUSTER_ID"             
#> [23] "MSLEVEL"                

## The variable "TITLE" represents the title of the spectrum defined in the
## MGF file
be$TITLE
#>  [1] "File193 Spectrum1719 scans: 2162"                                                         
#>  [2] "File193 Spectrum1944 scans: 2406"                                                         
#>  [3] "File193 Spectrum1968 scans: 2432"                                                         
#>  [4] "mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840:WNQLQAFWGTGK/2"
#>  [5] "mzspec:PXD002084:TCGA-AA-A01D-01A-23_W_VU_20121106_A0218_5I_R_FR15:scan:5209:DLTDYLMK/2"  
#>  [6] "mzspec:MSV000080679:j11962_C1orf144:scan:10671:DLTDYLMK/2"                                
#>  [7] "CCMSLIB00000840351"                                                                       
#>  [8] "blank_2-A,1_01_29559.812.812.1"                                                           
#>  [9] "File193 Spectrum1719 scans: 2162"                                                         
#> [10] "File193 Spectrum1719 scans: 2162"                                                         

## The default mapping of MGF fields to spectra variables is provided by
## the spectraVariableMapping function
spectraVariableMapping(MsBackendMgf())
#>              rtime     acquisitionNum        precursorMz precursorIntensity 
#>      "RTINSECONDS"            "SCANS"          "PEPMASS"       "PEPMASSINT" 
#>    precursorCharge 
#>           "CHARGE" 

## We can provide our own mapping e.g. to map the MGF field "TITLE" to a
## variable named "spectrumName":
map <- c(spectrumName = "TITLE", spectraVariableMapping(MsBackendMgf()))
map
#>       spectrumName              rtime     acquisitionNum        precursorMz 
#>            "TITLE"      "RTINSECONDS"            "SCANS"          "PEPMASS" 
#> precursorIntensity    precursorCharge 
#>       "PEPMASSINT"           "CHARGE" 

## We can then pass this mapping with parameter `mapping` to the
## backendInitialize method:
be <- backendInitialize(MsBackendMgf(), fls, mapping = map)
#> Start data import from 4 files ... 
#> done

## The title is now available as variable named spectrumName
be$spectrumName
#>  [1] "File193 Spectrum1719 scans: 2162"                                                         
#>  [2] "File193 Spectrum1944 scans: 2406"                                                         
#>  [3] "File193 Spectrum1968 scans: 2432"                                                         
#>  [4] "mzspec:PXD004732:01650b_BC2-TUM_first_pool_53_01_01-3xHCD-1h-R2:scan:41840:WNQLQAFWGTGK/2"
#>  [5] "mzspec:PXD002084:TCGA-AA-A01D-01A-23_W_VU_20121106_A0218_5I_R_FR15:scan:5209:DLTDYLMK/2"  
#>  [6] "mzspec:MSV000080679:j11962_C1orf144:scan:10671:DLTDYLMK/2"                                
#>  [7] "CCMSLIB00000840351"                                                                       
#>  [8] "blank_2-A,1_01_29559.812.812.1"                                                           
#>  [9] "File193 Spectrum1719 scans: 2162"                                                         
#> [10] "File193 Spectrum1719 scans: 2162"                                                         

## Next we create a Spectra object with this data
sps <- Spectra(be)

## We can use the 'MsBackendMgf' also to export spectra data in mgf format.
out_file <- tempfile()
export(sps, backend = MsBackendMgf(), file = out_file, map = map)

## The first 20 lines of the generated file:
readLines(out_file, n = 20)
#>  [1] "BEGIN IONS"                            
#>  [2] "msLevel=2"                             
#>  [3] "RTINSECONDS=1028"                      
#>  [4] "SCANS=2162"                            
#>  [5] "centroided=TRUE"                       
#>  [6] "PEPMASS=816.33826"                     
#>  [7] "CHARGE=2+"                             
#>  [8] "TITLE=File193 Spectrum1719 scans: 2162"
#>  [9] "102.0548 753.738"                      
#> [10] "103.00494 385.376"                     
#> [11] "103.03531 315.441"                     
#> [12] "115.05001 413.206"                     
#> [13] "115.08686 588.273"                     
#> [14] "120.08063 800.016"                     
#> [15] "124.10555 526.761"                     
#> [16] "124.20652 326.013"                     
#> [17] "1369.57385 1485.17"                    
#> [18] "1370.57495 1280.2"                     
#> [19] "1387.59021 3038.73"                    
#> [20] "1388.58691 2016.43"                    

## Next we add a new spectra variable to each spectrum
sps$spectrum_idx <- seq_along(sps)

## This new spectra variable will also be exported to the mgf file:
export(sps, backend = MsBackendMgf(), file = out_file, map = map)
readLines(out_file, n = 20)
#>  [1] "BEGIN IONS"                            
#>  [2] "msLevel=2"                             
#>  [3] "RTINSECONDS=1028"                      
#>  [4] "SCANS=2162"                            
#>  [5] "centroided=TRUE"                       
#>  [6] "PEPMASS=816.33826"                     
#>  [7] "CHARGE=2+"                             
#>  [8] "TITLE=File193 Spectrum1719 scans: 2162"
#>  [9] "spectrum_idx=1"                        
#> [10] "102.0548 753.738"                      
#> [11] "103.00494 385.376"                     
#> [12] "103.03531 315.441"                     
#> [13] "115.05001 413.206"                     
#> [14] "115.08686 588.273"                     
#> [15] "120.08063 800.016"                     
#> [16] "124.10555 526.761"                     
#> [17] "124.20652 326.013"                     
#> [18] "1369.57385 1485.17"                    
#> [19] "1370.57495 1280.2"                     
#> [20] "1387.59021 3038.73"                    

####
## Annotated MGF

## An example of a supported annotated MGF file
fl <- system.file("extdata", "xfiora.mgf", package = "MsBackendMgf")

## Lines with peak data start with a numeric and information is
## separated by a whitespace. The first two elements are the peak's m/z
## and intensity while any additional information is considered as
## annotation. Information for each peak is provided in one line.
readLines(fl)
#>   [1] "BEGIN IONS"                                                                                                                                 
#>   [2] "TITLE=PC(PGF1alpha/22:6(4Z,7Z,10Z,13Z,16Z,19Z))"                                                                                            
#>   [3] "SMILES=CCCCCC@@HO"                                                                                                                          
#>   [4] "PRECURSORTYPE=[M+H]+"                                                                                                                       
#>   [5] "COLLISIONENERGY=10.0"                                                                                                                       
#>   [6] "INSTRUMENTTYPE=HCD"                                                                                                                         
#>   [7] "COMMENT=\"In silico generated spectrum by Fiora OS v0.1.0\""                                                                                
#>   [8] "906.5854755900899 0.24612045288085938 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                 
#>   [9] "57.06987670808999 0.0012485487386584282 CCCC//[M-H]+"                                                                                       
#>  [10] "71.08552677208999 0.0015523125184699893 CCCCC//[M-H]+"                                                                                      
#>  [11] "100.08826642409 0.0031058108434081078 CCCCC[CH]O//[M-H]+"                                                                                   
#>  [12] "804.48101052209 0.003057875670492649 C=C[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                                    
#>  [13] "888.575459486 0.002360690850764513 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M]+"                      
#>  [14] "887.56708587409 0.02812790870666504 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                   
#>  [15] "790.4653604580899 0.0010008384706452489 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"              
#>  [16] "887.56708587409 0.002026958856731653 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                  
#>  [17] "19.01784113609 0.009557841811329126 O//[M+H]+"                                                                                              
#>  [18] "887.56708587409 0.0032752526458352804 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                 
#>  [19] "678.4129309620901 0.0013166849967092276 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#>  [20] "311.25807140009 0.0020137864630669355 CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M-H]+"                                                                
#>  [21] "339.25298602008996 0.03245394676923752 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                        
#>  [22] "568.33976602209 0.14936257898807526 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                   
#>  [23] "357.26355070408994 0.003117065876722336 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M+H]+"                                                       
#>  [24] "550.32920133809 0.012815228663384914 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                  
#>  [25] "578.3457939059999 0.002960008103400469 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M]+"                                                          
#>  [26] "577.33742029409 0.030111582949757576 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                          
#>  [27] "329.24750671609 0.002531866542994976 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)O//[M+H]+"                                           
#>  [28] "725.5350810960899 0.0012065907940268517 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M+H]+"
#>  [29] "724.527804644 0.0018646189710125327 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M]+"      
#>  [30] "723.51943103209 0.038313210010528564 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M-H]+"   
#>  [31] "183.06604455800002 0.005070459563285112 CN+(C)CCOP(=O)([O-])O//[M]+"                                                                        
#>  [32] "15.995463199909999 0.004359378479421139 [O-]//[M]+"                                                                                         
#>  [33] "104.10699048808999 0.0013647250598296523 CN+(C)CCO//[M]+"                                                                                   
#>  [34] "819.4817733339099 0.001787053421139717 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])O//[M]+"                          
#>  [35] "88.11207586808999 0.00103498506359756 CCN+(C)C//[M]+"                                                                                       
#>  [36] "87.10370225617999 0.008224652148783207 CCN+(C)C//[M-H]+"                                                                                    
#>  [37] "86.09642580408999 0.011802508495748043 CCN+(C)C//[M-2H]+"                                                                                   
#>  [38] "85.08805219217999 0.0019315838580951095 CCN+(C)C//[M-3H]+"                                                                                  
#>  [39] "73.08805219217999 0.002003903966397047 CN+(C)C//[M-H]+"                                                                                     
#>  [40] "846.50469985 0.003967077471315861 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCC//[M-H]+"                           
#>  [41] "60.08022716018 0.002686966210603714 CN+C//[M+H]+"                                                                                           
#>  [42] "58.06457709618 0.0016727085458114743 CN+C//[M-H]+"                                                                                          
#>  [43] "596.35581001009 0.055958326905965805 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M+H]+"                                                          
#>  [44] "311.23694203209004 0.015817463397979736 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC=O//[M-H]+"                                           
#>  [45] "283.24202741209 0.0027298536151647568 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                
#>  [46] "255.21072728408998 0.0011316784657537937 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                  
#>  [47] "662.3663746940899 0.0015134515706449747 CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-3H]+"                                                      
#>  [48] "241.19507722008998 0.0017022188985720277 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                   
#>  [49] "676.38202475809 0.0012113102711737156 C=CCCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                         
#>  [50] "674.3663746940899 0.0030237638857215643 C=CCCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                      
#>  [51] "229.19507722008998 0.004923875909298658 C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                      
#>  [52] "690.39767482209 0.007068068254739046 C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                      
#>  [53] "688.38202475809 0.0011182071175426245 C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                    
#>  [54] "217.19507722008998 0.004412529990077019 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M+H]+"                                                          
#>  [55] "215.17942715609 0.00425428943708539 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                              
#>  [56] "213.16377709209 0.0066286781802773476 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                           
#>  [57] "702.3976748220899 0.005614385008811951 CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                  
#>  [58] "201.16377709209 0.004419034346938133 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                              
#>  [59] "718.4289749500899 0.002718859352171421 C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                                 
#>  [60] "716.41332488609 0.002932032337412238 C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                   
#>  [61] "714.3976748220899 0.00470602186396718 C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                 
#>  [62] "189.16377709209 0.006993195042014122 C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                                
#>  [63] "187.14812702809 0.0015772052574902773 C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                              
#>  [64] "730.42897495009 0.006089580710977316 C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                               
#>  [65] "728.41332488609 0.0011209301883354783 C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                             
#>  [66] "177.16377709209 0.0047977156937122345 C=CC/C=C\\C/C=C\\C/C=C\\CC//[M+H]+"                                                                   
#>  [67] "175.14812702809 0.0039370497688651085 C=CC/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                                   
#>  [68] "173.13247696409 0.006436163559556007 C=CC/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                   
#>  [69] "742.4289749500899 0.0065302494913339615 CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                          
#>  [70] "161.13247696409 0.005028888117522001 CC/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                     
#>  [71] "758.4602750780899 0.0027850852347910404 C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                         
#>  [72] "756.44462501409 0.0030861366540193558 C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                           
#>  [73] "754.4289749500899 0.004762527532875538 C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                         
#>  [74] "149.13247696409 0.007314636372029781 C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                                       
#>  [75] "147.11682690009 0.0016508938279002905 C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                     
#>  [76] "770.46027507809 0.006089580710977316 C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                        
#>  [77] "768.44462501409 0.0011209301883354783 C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                      
#>  [78] "137.13247696409 0.0047977156937122345 C=CC/C=C\\C/C=C\\CC//[M+H]+"                                                                          
#>  [79] "135.11682690009 0.0039370497688651085 C=CC/C=C\\C/C=C\\CC//[M-H]+"                                                                          
#>  [80] "133.10117683609 0.006436163559556007 C=CC/C=C\\C/C=C\\CC//[M-3H]+"                                                                          
#>  [81] "782.4602750780899 0.0065302494913339615 CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                   
#>  [82] "121.10117683608999 0.005028888117522001 CC/C=C\\C/C=C\\CC//[M-3H]+"                                                                         
#>  [83] "798.4915752060899 0.0027850852347910404 C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                  
#>  [84] "796.47592514209 0.0030861366540193558 C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                    
#>  [85] "794.4602750780899 0.004762527532875538 C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                  
#>  [86] "109.10117683608999 0.007314636372029781 C/C=C\\C/C=C\\CC//[M-H]+"                                                                           
#>  [87] "107.08552677208999 0.0016508938279002905 C/C=C\\C/C=C\\CC//[M-3H]+"                                                                         
#>  [88] "810.49157520609 0.006118484307080507 C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                 
#>  [89] "808.47592514209 0.0010935019236057997 C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                               
#>  [90] "97.10117683608999 0.00479193776845932 C=CC/C=C\\CC//[M+H]+"                                                                                 
#>  [91] "95.08552677208999 0.00393778458237648 C=CC/C=C\\CC//[M-H]+"                                                                                 
#>  [92] "93.06987670808999 0.00646106107160449 C=CC/C=C\\CC//[M-3H]+"                                                                                
#>  [93] "822.4915752060899 0.007105632219463587 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                             
#>  [94] "81.06987670808999 0.0054127830080688 CC/C=C\\CC//[M-3H]+"                                                                                   
#>  [95] "838.5228753340899 0.0024260678328573704 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                           
#>  [96] "836.50722527009 0.003049242775887251 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                              
#>  [97] "834.4915752060899 0.004425586201250553 C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                           
#>  [98] "69.06987670808999 0.008552799001336098 C/C=C\\CC//[M-H]+"                                                                                   
#>  [99] "67.05422664409 0.0017731321277096868 C/C=C\\CC//[M-3H]+"                                                                                    
#> [100] "850.52287533409 0.00296195806004107 C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                           
#> [101] "57.06987670808999 0.002908135997131467 C=CCC//[M+H]+"                                                                                       
#> [102] "55.054226644089994 0.0025195078924298286 C=CCC//[M-H]+"                                                                                     
#> [103] "53.038576580089995 0.0017213564133271575 C=CCC//[M-3H]+"                                                                                    
#> [104] "862.5228753340899 0.0011330428533256054 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                     
#> [105] "END IONS"                                                                                                                                   

## Importing the data using an `MsBackendAnnotatedMgf`
ba <- backendInitialize(MsBackendAnnotatedMgf(), fl)
#> Start data import from 1 files ... 
#> done
ba
#> MsBackendAnnotatedMgf with 1 spectra
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         2        NA        NA
#>  ... 23 more variables/columns.

## An additional peaks variable is available.
peaksVariables(ba)
#> [1] "mz"        "intensity" "V1"       

ba$V1
#> [[1]]
#>  [1] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"               
#>  [2] "CCCC//[M-H]+"                                                                                       
#>  [3] "CCCCC//[M-H]+"                                                                                      
#>  [4] "CCCCC[CH]O//[M-H]+"                                                                                 
#>  [5] "C=C[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                                 
#>  [6] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M]+"                 
#>  [7] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#>  [8] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"              
#>  [9] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#> [10] "O//[M+H]+"                                                                                          
#> [11] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#> [12] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#> [13] "CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M-H]+"                                                              
#> [14] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                       
#> [15] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"               
#> [16] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M+H]+"                                                       
#> [17] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"               
#> [18] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M]+"                                                         
#> [19] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+"                                                       
#> [20] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)O//[M+H]+"                                        
#> [21] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M+H]+"
#> [22] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M]+"  
#> [23] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@@HCOC(=O)CCCCCC[C@@h]1C@@HC@HC[C@@h]1O//[M-H]+"
#> [24] "CN+(C)CCOP(=O)([O-])O//[M]+"                                                                        
#> [25] "[O-]//[M]+"                                                                                         
#> [26] "CN+(C)CCO//[M]+"                                                                                    
#> [27] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])O//[M]+"                         
#> [28] "CCN+(C)C//[M]+"                                                                                     
#> [29] "CCN+(C)C//[M-H]+"                                                                                   
#> [30] "CCN+(C)C//[M-2H]+"                                                                                  
#> [31] "CCN+(C)C//[M-3H]+"                                                                                  
#> [32] "CN+(C)C//[M-H]+"                                                                                    
#> [33] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCC//[M-H]+"                     
#> [34] "CN+C//[M+H]+"                                                                                       
#> [35] "CN+C//[M-H]+"                                                                                       
#> [36] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M+H]+"                                                       
#> [37] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC=O//[M-H]+"                                           
#> [38] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                              
#> [39] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                   
#> [40] "CCCCCC@H/C=C/[C@@h]1C@@HC@@HC[C@H]1O//[M-3H]+"                                                      
#> [41] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                    
#> [42] "C=CCCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                       
#> [43] "C=CCCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                      
#> [44] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                      
#> [45] "C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                   
#> [46] "C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                  
#> [47] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M+H]+"                                                          
#> [48] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                          
#> [49] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                         
#> [50] "CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                                 
#> [51] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                           
#> [52] "C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                                
#> [53] "C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                                
#> [54] "C=CC/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                               
#> [55] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                             
#> [56] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                            
#> [57] "C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                            
#> [58] "C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                           
#> [59] "C=CC/C=C\\C/C=C\\C/C=C\\CC//[M+H]+"                                                                 
#> [60] "C=CC/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                                 
#> [61] "C=CC/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                
#> [62] "CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                          
#> [63] "CC/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                  
#> [64] "C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                         
#> [65] "C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                         
#> [66] "C=CC/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                        
#> [67] "C/C=C\\C/C=C\\C/C=C\\CC//[M-H]+"                                                                    
#> [68] "C/C=C\\C/C=C\\C/C=C\\CC//[M-3H]+"                                                                   
#> [69] "C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                     
#> [70] "C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                    
#> [71] "C=CC/C=C\\C/C=C\\CC//[M+H]+"                                                                        
#> [72] "C=CC/C=C\\C/C=C\\CC//[M-H]+"                                                                        
#> [73] "C=CC/C=C\\C/C=C\\CC//[M-3H]+"                                                                       
#> [74] "CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                   
#> [75] "CC/C=C\\C/C=C\\CC//[M-3H]+"                                                                         
#> [76] "C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                                  
#> [77] "C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                                  
#> [78] "C=CC/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                                 
#> [79] "C/C=C\\C/C=C\\CC//[M-H]+"                                                                           
#> [80] "C/C=C\\C/C=C\\CC//[M-3H]+"                                                                          
#> [81] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                              
#> [82] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                             
#> [83] "C=CC/C=C\\CC//[M+H]+"                                                                               
#> [84] "C=CC/C=C\\CC//[M-H]+"                                                                               
#> [85] "C=CC/C=C\\CC//[M-3H]+"                                                                              
#> [86] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                            
#> [87] "CC/C=C\\CC//[M-3H]+"                                                                                
#> [88] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+"                           
#> [89] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                           
#> [90] "C=CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                          
#> [91] "C/C=C\\CC//[M-H]+"                                                                                  
#> [92] "C/C=C\\CC//[M-3H]+"                                                                                 
#> [93] "C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-H]+"                       
#> [94] "C=CCC//[M+H]+"                                                                                      
#> [95] "C=CCC//[M-H]+"                                                                                      
#> [96] "C=CCC//[M-3H]+"                                                                                     
#> [97] "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M-3H]+"                     
#> 

## The length of such peaks variables is the same as the length of the
## m/z or intensity values, i.e. each peak has one value (with the value
## being `NA` if missing).
length(ba$V1[[1L]])
#> [1] 97
length(ba$mz[[1L]])
#> [1] 97

## Extracting the peaks data from a `Spectra` with a `MsBackendAnnotatedMgf`
s <- Spectra(ba)
pd <- peaksData(s, peaksVariables(ba))[[1L]]
head(pd)
#>         mz   intensity
#> 1 15.99546 0.004359378
#> 2 19.01784 0.009557842
#> 3 53.03858 0.001721356
#> 4 55.05423 0.002519508
#> 5 57.06988 0.001248549
#> 6 57.06988 0.002908136
#>                                                                                     V1
#> 1 CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M+H]+
#> 2                                                                         CCCC//[M-H]+
#> 3                                                                        CCCCC//[M-H]+
#> 4                                                                   CCCCC[CH]O//[M-H]+
#> 5                                                   C=C[C@@h]1C@@HC@@HC[C@H]1O//[M-H]+
#> 6   CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)OC@HCOP(=O)([O-])OCCN+(C)C//[M]+
class(pd)
#> [1] "data.frame"
```
