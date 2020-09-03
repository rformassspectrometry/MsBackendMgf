##' @param f `character(1)` with the path to an mgf file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param mapping named `character` vector to rename mgf fields to spectra
##'     variables.
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame
##'
##' @importFrom IRanges NumericList
##'
##' @importFrom MsCoreUtils rbindFill
##'
##' @author Laurent Gatto, Johannes Rainer
##'
##' @noRd
.read_mgf <- function(f, msLevel = 2L,
                      mapping = spectraVariableMapping(), ...) {
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
                                         mapping = mapping)

    res <- DataFrame(rbindFill(sp))

    spv <- Spectra:::.SPECTRA_DATA_COLUMNS
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- IRanges::NumericList(res$mz)
    res$intensity <- IRanges::NumericList(res$intensity)
    res$dataOrigin <- f
    res$msLevel <- as.integer(msLevel)
    res
}

##' @description
##'
##' Extract **all** fields from the MGF eventually renaming the field names to
##' the spectra variable names specified with `mapping`.
##'
##' @param mgf `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @param mapping named `character` providing the mapping of MGF fields to
##'     spectra variable names.
##'
##' @author Laurent Gatto, Johannes Rainer
##'
##' @importFrom stats setNames
##'
##' @noRd
.extract_mgf_spectrum <- function(mgf, mapping = spectraVariableMapping()) {
    ## grep description
    desc.idx <- grep("=", mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]

    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
    mode(ms) <- "double"

    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)),
                     substring(desc, 1L, r - 1L))
    title <- unname(desc["TITLE"])

    desc[c("PEPMASS", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]

    ## Use all fields in the MGF renaming the ones specified by mapping.
    desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])

    idx <- match(names(desc), mapping)
    not_na <- !is.na(idx)
    if (any(not_na))
        names(desc)[not_na] <- names(mapping)[idx][not_na]
    res <- data.frame(matrix(desc, nrow = 1,
                             dimnames = list(NULL, names(desc))))
    res$mz = list(ms[, 1L])
    res$intensity = list(ms[, 2L])
    res
}
