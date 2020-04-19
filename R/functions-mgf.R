##' @importFrom S4Vectors DataFrame
##'
##' @importFrom IRanges NumericList
##' 
##' @author Laurent Gatto
##' 
##' @noRd
##' 
##' @param f `character(1)` holding the path to an mgf file.
##' 
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##' 
##' @param keyValues `data.frame` with mgf key/values
##'     pairs. Default is [mgfKeyValues()].
##' 
##' @param ... Additional parameters, currently ignored.
##'
.read_mgf <- function(f, msLevel = 2L,
                      keyValues = mgfKeyValues(),
                      ...) {
    if (length(f) != 1L)
        stop("Please provide a single mgf file.")

    mgf <- scan(file = f, what = "",
                sep = "\n", quote = "",
                allowEscapes = FALSE,
                quiet = TRUE)
    ## From http://www.matrixscience.com/help/data_file_help.html#GEN
    ## Comment lines beginning with one of the symbols #;!/ can be
    ## included, but only outside of the BEGIN IONS and END IONS
    ## statements that delimit an MS/MS dataset.
    cmts <- grep("^[#;!/]", mgf)
    if (length(cmts))
        mgf <- mgf[-cmts]

    begin <- grep("BEGIN IONS", mgf) + 1L
    end <- grep("END IONS", mgf) - 1L
    n <- length(begin)
    sp <- vector("list", length = n)

    for (i in seq(along = sp)) 
        sp[[i]] <- .extract_mgf_spectrum(mgf[begin[i]:end[i]],
                                         keyValues)
    res <- DataFrame(do.call(rbind, sp))

    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
    }

    if (!"msLevel" %in% names(res))
        res$msLevel <- as.integer(msLevel)

    res$mz <- IRanges::NumericList(res$mz)
    res$intensity <- IRanges::NumericList(res$intensity)
    res$precursorCharge <- as.integer(res$precursorCharge)
    res$scanIndex <- as.integer(res$scanIndex)
    res$dataOrigin <- f
    res
}


##' @param mgf `character()` containing lines from a single spectrum
##'     from an mgf file.
##' 
##' @param mgf_key_values `data.frame` containing key/value pairs to
##'     format the header keys.
##' 
##' @author Laurent Gatto
##' 
##' @importFrom stats setNames
##'
##' @noRd
.extract_mgf_spectrum <- function(mgf, mgf_key_values) {
   ## grep description
    desc.idx <- grep("=", mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]

    ## extract ms data
    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
    mode(ms) <- "double"

    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)

    ## extract header data
    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)),
                     substring(desc, 1L, r - 1L))    
    desc <- as.list(desc)

    browser()
    
    if ("PEPMASS" %in% names(desc)) {
        ## PEPMASSMZ and PEPMASSINT can contain 2 numericals,
        ## corresponding to the precursorMz and precursorIntensity
        pepmass <- as.numeric(strsplit(desc[["PEPMASS"]], "[[:space:]]+")[[1L]])
        if (length(pepmass) == 1) 
            pepmass <- c(pepmass, NA)
        ## assuming length of 2, additional values ignored
        desc[c("precursorMz", "precursorIntensity")] <- pepmass[1:2]
        desc[["PEPMASS"]] <- NULL 
    }        

    for (i in seq_along(desc)) {
        key <- names(desc)[i]
        if (key %in% mgf_key_values[[1]]) {
            names(desc)[i]  <- mgf_key_values[mgf_key_values[[1]] == key, 2]
            ## CHARGE sometimes ends with a '+' or '-'
            desc[[i]] <- sub("[+-]", "", desc[[i]])
            desc[[i]] <- as.numeric(desc[[i]])            
        }
    }
    
    c(desc[order(names(desc))],
      mz = list(ms[, 1L]),
      intensity = list(ms[, 2L]))

}

##' It is possible to 
##'
##' @title Mgf key/value pairs
##' 
##' @return A `data.frame` with mgf key/value pairs.
##' 
##' @author Laurent Gatto
##'
##' @importFrom utils read.csv
##'
##' @export
##'
##' @examples
##' mgfKeyValues()
mgfKeyValues <- function() 
    read.csv(dir(system.file("extdata",
                             package = "MsBackendMgf"),
                 pattern = "mgf_key_values.csv",
                 full.names = TRUE),
             header = TRUE,
             stringsAsFactors = FALSE)
