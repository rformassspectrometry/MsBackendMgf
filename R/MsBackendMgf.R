#' @include hidden_aliases.R
NULL

#' @title MS data backend for MGF files
#'
#' @aliases MsBackendMgf-class
#'
#' @description
#'
#' The `MsBackendMgf` class supports import and export of MS/MS spectra data
#' from/to files in Mascot Generic Format
#' ([mgf](http://www.matrixscience.com/help/data_file_help.html))
#' files. After initial import, the full MS data is kept in
#' memory. `MsBackendMgf` extends the [Spectra::MsBackendDataFrame()] backend
#' directly and supports thus the [Spectra::applyProcessing()] function to make
#' data manipulations persistent.
#'
#' The `MsBackendAnnotatedMgf` class supports import of data from MGF files
#' that provide, in addition to the *m/z* and intensity values, also
#' additional annotations/metadata for each mass peak. For such MGF files it
#' is expected that each line contains information from a single mass peak,
#' separated by a white space (blank). The first two elements are expected to
#' be the peak's *m/z* and intensity values, while each additional element is
#' considered an annotation for this specific peak. See examples below for the
#' format of a supported MGF file. The `backendInitialize()` method of
#' `MsBackendAnnotatedMgf` does not support parameter `nlines`. Also, import
#' of data can be considerably slower compared to the standard `MsBackendMgf`
#' backend, because of the additionally required parsing of peak annotations.
#' Peaks information in MGF files are not named, thus, additional peaks
#' annotations are named using the standard naming convention for column named
#' of data frames: the first peaks annotation is called `"V1"`, the second (if
#' available) `"V2"` and so on.
#'
#' New objects are created with the `MsBackendMgf()` or
#' `MsBackendAnnotatedMgf()` function. The `backendInitialize()` method has to
#' be subsequently called to initialize the object and import the MS/MS data
#' from (one or more) MGF files.
#'
#' The `MsBackendMgf` backend provides an `export` method that allows to export
#' the data from the `Spectra` object (parameter `x`) to a file in mgf format.
#' See the package vignette for details and examples.
#'
#' Default mappings from fields in the MGF file to spectra variable names are
#' provided by the `spectraVariableMapping()` function. This function returns a
#' named character vector were names are the spectra variable names and the
#' values the respective field names in the MGF files. This named character
#' vector is submitted to the import and export function with parameter
#' `mapping`. It is also possible to pass own mappings (e.g. for special
#' MGF dialects) with the `mapping` parameter.
#'
#' @param object Instance of `MsBackendMgf` class.
#'
#' @param file `character(1)` with the (full) file name to which the data
#'     should be exported.
#'
#' @param files `character` with the (full) file name(s) of the mgf file(s)
#'     from which MS/MS data should be imported.
#'
#' @param format for `spectraVariableMapping()`: `character(1)` defining the
#'     format to be used. Currently only `format = "mgf"` is supported.
#'
#' @param mapping for `backendInitialize()` and `export`: named `character`
#'     vector allowing to specify how fields from the MGF file should be
#'     renamed. Names are supposed to be the spectra variable name and
#'     values of the vector the field names in the MGF file. See output of
#'     `spectraVariableMapping()` for the expected format and examples
#'     below or description above for details.
#'
#' @param nlines for `backendInitialize()` of `MsBackendMgf`: `integer(1)`
#'     defining the number of lines that should be imported and processed
#'     from the MGF file(s).
#'     By default (`nlines = -1L`) the full file is imported and processed at
#'     once. If set to a positive integer, the data is imported and processed
#'     *chunk-wise* using [readMgfSplit()].
#'
#' @param exportTitle `logical(1)` whether the *TITLE* field should be included
#'     in the exported MGF file. If `TRUE` (the default) a `spectraVariable`
#'     called `"TITLE"` will be used, if no such variable is present either the
#'     `spectraNames(object)` will be used or, if they are empty, a title will
#'     be generated including the MS level, retention time and acquisition
#'     number of the spectrum.
#'
#' @param x for `export()`: an instance of [Spectra::Spectra()] class with the
#'     data that should be exported.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup. If parallel processing is enabled (with `BPPARAM` different than
#'     `SerialParam()`, the default) and length of `files` is larger than one,
#'     import is performed in parallel on a per-file basis. If data is to be
#'     imported from a single file (i.e., length of `files` is one), parsing
#'     of the imported file is performed in parallel. See also
#'     [BiocParallel::SerialParam()] for information on available parallel
#'     processing setup options.
#'
#' @param ... Currently ignored.
#'
#' @author Laurent Gatto, Corey Broeckling and Johannes Rainer
#'
#' @importClassesFrom Spectra MsBackendDataFrame
#'
#' @exportClass MsBackendMgf
#'
#' @name MsBackendMgf
#'
#' @return
#'
#' See description above.
#'
#' @examples
#'
#' library(BiocParallel)
#' fls <- dir(system.file("extdata", package = "MsBackendMgf"),
#'     full.names = TRUE, pattern = "mgf$")
#'
#' ## Create an MsBackendMgf backend and import data from test mgf files.
#' be <- backendInitialize(MsBackendMgf(), fls)
#' be
#'
#' be$msLevel
#' be$intensity
#' be$mz
#'
#' ## The spectra variables that are available; note that not all of them
#' ## have been imported from the MGF files.
#' spectraVariables(be)
#'
#' ## The variable "TITLE" represents the title of the spectrum defined in the
#' ## MGF file
#' be$TITLE
#'
#' ## The default mapping of MGF fields to spectra variables is provided by
#' ## the spectraVariableMapping function
#' spectraVariableMapping(MsBackendMgf())
#'
#' ## We can provide our own mapping e.g. to map the MGF field "TITLE" to a
#' ## variable named "spectrumName":
#' map <- c(spectrumName = "TITLE", spectraVariableMapping(MsBackendMgf()))
#' map
#'
#' ## We can then pass this mapping with parameter `mapping` to the
#' ## backendInitialize method:
#' be <- backendInitialize(MsBackendMgf(), fls, mapping = map)
#'
#' ## The title is now available as variable named spectrumName
#' be$spectrumName
#'
#' ## Next we create a Spectra object with this data
#' sps <- Spectra(be)
#'
#' ## We can use the 'MsBackendMgf' also to export spectra data in mgf format.
#' out_file <- tempfile()
#' export(sps, backend = MsBackendMgf(), file = out_file, map = map)
#'
#' ## The first 20 lines of the generated file:
#' readLines(out_file, n = 20)
#'
#' ## Next we add a new spectra variable to each spectrum
#' sps$spectrum_idx <- seq_along(sps)
#'
#' ## This new spectra variable will also be exported to the mgf file:
#' export(sps, backend = MsBackendMgf(), file = out_file, map = map)
#' readLines(out_file, n = 20)
#'
#' ####
#' ## Annotated MGF
#'
#' ## An example of a supported annotated MGF file
#' fl <- system.file("extdata", "xfiora.mgf", package = "MsBackendMgf")
#'
#' ## Lines with peak data start with a numeric and information is
#' ## separated by a whitespace. The first two elements are the peak's m/z
#' ## and intensity while any additional information is considered as
#' ## annotation. Information for each peak is provided in one line.
#' readLines(fl)
#'
#' ## Importing the data using an `MsBackendAnnotatedMgf`
#' ba <- backendInitialize(MsBackendAnnotatedMgf(), fl)
#' ba
#'
#' ## An additional peaks variable is available.
#' peaksVariables(ba)
#'
#' ba$V1
#'
#' ## The length of such peaks variables is the same as the length of the
#' ## m/z or intensity values, i.e. each peak has one value (with the value
#' ## being `NA` if missing).
#' length(ba$V1[[1L]])
#' length(ba$mz[[1L]])
#'
#' ## Extracting the peaks data from a `Spectra` with a `MsBackendAnnotatedMgf`
#' s <- Spectra(ba)
#' pd <- peaksData(s, peaksVariables(ba))[[1L]]
#' head(pd)
#' class(pd)
NULL

setClass("MsBackendMgf",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra spectraData<- $<- $
#'
#' @importFrom BiocParallel SerialParam
#'
#' @importMethodsFrom ProtGenerics backendInitialize
#'
#' @importMethodsFrom Spectra backendInitialize
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @importFrom methods validObject
#'
#' @exportMethod backendInitialize
#'
#' @rdname MsBackendMgf
setMethod("backendInitialize", signature = "MsBackendMgf",
          function(object, files, mapping = spectraVariableMapping(object),
                   nlines = -1L, ..., BPPARAM = SerialParam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for ", class(object))
              if (!is.character(files))
                  stop("Parameter 'files' is expected to be a character vector",
                       " with the files names from where data should be",
                       " imported")
              if (!is.numeric(nlines))
                  stop("'nlines' needs to be an integer")
              nlines <- as.integer(nlines)
              files <- normalizePath(files)
              if (any(!file.exists(files)))
                  stop("file(s) ",
                       paste(files[!file.exists(files)], collapse = ", "),
                       " not found")
              ## Import data and rbind.
              message("Start data import from ", length(files), " files ... ",
                      appendLF = FALSE)
              if (nlines > 0)
                  FUN <- readMgfSplit
              else FUN <- readMgf
              if (length(files) > 1) {
                  res <- bplapply(files, FUN = FUN, mapping = mapping,
                                  nlines = nlines, BPPARAM = BPPARAM)
                  res <- do.call(rbindFill, res)
              } else
                  res <- FUN(files, mapping = mapping, nlines = nlines,
                             BPPARAM = BPPARAM)
              message("done")
              res@metadata <- list()
              spectraData(object) <- res
              object$dataStorage <- "<memory>"
              object$centroided <- TRUE
              validObject(object)
              object
          })

#' @rdname MsBackendMgf
#'
#' @importFrom methods new
#'
#' @export MsBackendMgf
MsBackendMgf <- function() {
    new("MsBackendMgf")
}

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @exportMethod spectraVariableMapping
#'
#' @rdname MsBackendMgf
setMethod("spectraVariableMapping", "MsBackendMgf",
          function(object, format = c("mgf")) {
              switch(match.arg(format),
                     "mgf" = c(
                         rtime = "RTINSECONDS",
                         acquisitionNum = "SCANS",
                         precursorMz = "PEPMASS",
                         precursorIntensity = "PEPMASSINT",
                         precursorCharge = "CHARGE"
                     )
                     )
          })

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMgf
setMethod("export", "MsBackendMgf",
          function(object, x, file = tempfile(),
                   mapping = spectraVariableMapping(object),
                   exportTitle = TRUE, ...) {
              if (missing(x))
                  stop("Required parameter 'x' is missing. 'x' should be a ",
                       "'Spectra' object with the full spectra data.")
              if (!inherits(x, "Spectra"))
                  stop("Parameter 'x' is supposed to be a 'Spectra' object ",
                       "with the full spectra data to be exported.")
              .export_mgf(x = x, con = file, mapping = mapping,
                          exportTitle = exportTitle)
          })

setClass("MsBackendAnnotatedMgf",
         contains = "MsBackendMgf",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @rdname MsBackendMgf
#'
#' @importMethodsFrom S4Vectors metadata
setMethod("backendInitialize", signature = "MsBackendAnnotatedMgf",
          function(object, files, mapping = spectraVariableMapping(object),
                   ..., BPPARAM = SerialParam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for ", class(object))
              if (!is.character(files))
                  stop("Parameter 'files' is expected to be a character vector",
                       " with the files names from where data should be",
                       " imported")
              files <- normalizePath(files)
              if (any(!file.exists(files)))
                  stop("file(s) ",
                       paste(files[!file.exists(files)], collapse = ", "),
                       " not found")
              ## Import data and rbind.
              message("Start data import from ", length(files), " files ... ",
                      appendLF = FALSE)
              if (length(files) > 1) {
                  res <- bplapply(files, FUN = readMgf, mapping = mapping,
                                  annotated = TRUE, BPPARAM = BPPARAM)
                  pcol <- unique(unlist(lapply(res, metadata)))
                  res <- do.call(rbindFill, res)
              } else {
                  res <- readMgf(files, mapping = mapping, annotated = TRUE,
                                 BPPARAM = BPPARAM)
                  pcol <- unique(unlist(metadata(res)))
              }
              message("done")
              res@metadata <- list()
              spectraData(object) <- res
              object@peaksVariables <- c("mz", "intensity", pcol)
              object$dataStorage <- "<memory>"
              object$centroided <- TRUE
              validObject(object)
              object
          })

#' @rdname MsBackendMgf
#'
#' @export MsBackendAnnotatedMgf
MsBackendAnnotatedMgf <- function() {
    new("MsBackendAnnotatedMgf")
}
