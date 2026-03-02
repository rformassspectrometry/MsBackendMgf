# Reading MGF files

The `readMgf()` function imports the data from a file in MGF format
reading all specified fields and returning the data as a
[`S4Vectors::DataFrame()`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

For very large MGF files the `readMgfSplit()` function might be used
instead. In contrast to the `readMgf()` functions, `readMgfSplit()`
reads only `nlines` lines from an MGF file at once reducing thus the
memory demand (at the cost of a lower performance, compared to
`readMgf()`).

## Usage

``` r
readMgf(
  f,
  msLevel = 2L,
  mapping = spectraVariableMapping(MsBackendMgf()),
  annotated = FALSE,
  ...,
  BPPARAM = SerialParam()
)

readMgfSplit(
  f,
  msLevel = 2L,
  mapping = spectraVariableMapping(MsBackendMgf()),
  nlines = 1e+05,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- f:

  `character(1)` with the path to an mgf file.

- msLevel:

  `numeric(1)` with the MS level. Default is 2.

- mapping:

  named `character` vector to rename mgf fields to spectra variables.

- annotated:

  For `readMgf()`: `logical(1)` whether the MGF file contains additional
  peak annotations. See examples below or the documentation for
  [`MsBackendAnnotatedMgf()`](https://rformassspectrometry.github.io/MsBackendMgf/reference/MsBackendMgf.md)
  for information on the expected format.

- ...:

  Additional parameters, currently ignored.

- BPPARAM:

  parallel processing setup that should be used. Only the parsing of the
  imported MGF file is performed in parallel.

- nlines:

  for `readMgfSplit()`: `integer(1)` with the number of lines that
  should be imported and parsed in each iteration.

## Value

A `DataFrame` with each row containing the data from one spectrum in the
MGF file. m/z and intensity values are available in columns `"mz"` and
`"intensity"` in a list representation. For `readMgf()` with
`annotated = TRUE` also all peaks annotation columns (named
`"V1", etc) are provided in a list representation, with the lengths of elements matching those of `"mz"`or`"intensity"\`.

## Author

Laurent Gatto, Johannes Rainer, Sebastian Gibb, Corey Broeckling

## Examples

``` r

fls <- dir(system.file("extdata", package = "MsBackendMgf"),
    full.names = TRUE, pattern = "mgf$")[1L]

readMgf(fls)
#> DataFrame with 3 rows and 10 columns
#>                    TITLE precursorMz precursorCharge     rtime acquisitionNum
#>              <character>   <numeric>       <integer> <numeric>      <integer>
#> 1 File193 Spectrum1719..     816.338               2      1028           2162
#> 2 File193 Spectrum1944..     787.828               2      1117           2406
#> 3 File193 Spectrum1968..     780.840               2      1127           2432
#>   precursorIntensity                          mz                   intensity
#>            <numeric>               <NumericList>               <NumericList>
#> 1                 NA 102.055,103.005,103.035,... 753.738,385.376,315.441,...
#> 2             880650 101.071,102.055,103.002,... 1228.93,1424.66,1550.90,...
#> 3            1265631 102.056,103.000,115.051,... 1340.44,1714.76,1938.82,...
#>               dataOrigin   msLevel
#>              <character> <integer>
#> 1 /__w/_temp/Library/M..         2
#> 2 /__w/_temp/Library/M..         2
#> 3 /__w/_temp/Library/M..         2

## Annotated MGF
fl <- system.file("extdata", "xfiora.mgf", package = "MsBackendMgf")
res <- readMgf(fl, annotated = TRUE)
colnames(res)
#>  [1] "TITLE"              "SMILES"             "PRECURSORTYPE"     
#>  [4] "COLLISIONENERGY"    "INSTRUMENTTYPE"     "COMMENT"           
#>  [7] "precursorMz"        "precursorIntensity" "mz"                
#> [10] "intensity"          "V1"                 "dataOrigin"        
#> [13] "msLevel"           
res$V1
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
```
