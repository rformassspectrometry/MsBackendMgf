#' @include hidden_aliases.R
NULL

#' @title MS data backend for mgf files
#'
#' @aliases MsBackendMgf-class
#'
#' @description
#'
#' The `MsBackendMgf` class supports import and export of MS/MS spectra data
#' from/to files in Mascot Generic Format
#' ([mgf](http://www.matrixscience.com/help/data_file_help.html))
#' files. After initial import, the full MS data is kept in
#' memory. `MsBackendMgf` extends the [MsBackendDataFrame()] backend
#' directly and supports thus the [applyProcessing()] function to make
#' data manipulations persistent.
#'
#' New objects are created with the `MsBackendMgf` function. The
#' `backendInitialize` method has to be subsequently called to
#' initialize the object and import MS/MS data from (one or more) mgf
#' files.
#'
#' The `MsBackendMgf` backend provides an `export` method that allows to export
#' the data from the `Spectra` object (parameter `x`) to a file in mgf format.
#' See the package vignette for details and examples.
#'
#' Default mappings from fields in the MGF file to spectra variable names are
#' provided by the `spectraVariableMapping` function. This function returns a
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
#' @param format for `spectraVariableMapping`: `character(1)` defining the
#'     format to be used. Currently only `format = "mgf"` is supported.
#'
#' @param mapping for `backendInitialize` and `export`: named `character` vector
#'     allowing to specify how fields from the MGF file should be renamed. Names
#'     are supposed to be the spectra variable name and values of the vector
#'     the field names in the MGF file. See output of `spectraVariableMapping()`
#'     for the expected format and examples below or description above for
#'     details.
#'
#' @param x for `export`: an instance of [Spectra()] class with the data that
#'     should be exported.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     bpparam()`. See [bpparam()] for more information.
#'
#' @param ... Currently ignored.
#'
#' @author Laurent Gatto and Johannes Rainer
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
#' ## Create an MsBackendMgf backend and import data from test mgf files.
#' fls <- dir(system.file("extdata", package = "MsBackendMgf"),
#'     full.names = TRUE, pattern = "mgf$")
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
#' spectraVariableMapping()
#'
#' ## We can provide our own mapping e.g. to map the MGF field "TITLE" to a
#' ## variable named "spectrumName":
#' map <- c(spectrumName = "TITLE", spectraVariableMapping())
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
NULL

setClass("MsBackendMgf",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra backendInitialize spectraData<- $<- $
#'
#' @importFrom BiocParallel bpparam
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @importFrom methods validObject
#'
#' @exportMethod backendInitialize
#'
#' @rdname MsBackendMgf
setMethod("backendInitialize", signature = "MsBackendMgf",
          function(object, files, mapping = spectraVariableMapping(),
                   ..., BPPARAM = bpparam()) {
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
              res <- bplapply(files, FUN = .read_mgf,
                              mapping = mapping,
                              BPPARAM = BPPARAM)
              message("done")
              res <- do.call(rbindFill, res)
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

#' @export
#'
#' @rdname MsBackendMgf
spectraVariableMapping <- function(format = c("mgf")) {
    ## In future eventually define that in a text file and import upon package
    ## init.
    switch(match.arg(format),
           "mgf" = c(
               rtime = "RTINSECONDS",
               acquisitionNum = "SCANS",
               precursorMz = "PEPMASS",
               precursorIntensity = "PEPMASSINT",
               precursorCharge = "CHARGE"
           )
           )
}

#' @importMethodsFrom Spectra export
#'
#' @exportMethod export
#'
#' @rdname MsBackendMgf
setMethod("export", "MsBackendMgf", function(object, x, file = tempfile(),
                                             mapping = spectraVariableMapping(),
                                             ...) {
    if (missing(x))
        stop("Required parameter 'x' is missing. 'x' should be a 'Spectra' ",
             "object with the full spectra data.")
    if (!inherits(x, "Spectra"))
        stop("Parameter 'x' is supposed to be a 'Spectra' object with the full",
             " spectra data to be exported.")
    .export_mgf(x = x, con = file, mapping = mapping)
})
