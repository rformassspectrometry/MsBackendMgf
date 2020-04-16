##' @param f
##' 
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##' 
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame
##' 
##' @author Laurent Gatto
##' 
##' @noRd
.read_mgf <- function(f, msLevel = 2L, ...) {
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
        sp[[i]] <- .extract_mgf_spectrum(mgf[begin[i]:end[i]])

    res <- DataFrame(do.call(rbind, sp))

    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
    }

    res$dataOrigin <- f
    res$msLevel <- as.integer(msLevel)
    res
}

##' @param mgf `character()` of lines defining a spectrum in mgf
##'     format.
##' 
##' @author Laurent Gatto
##' 
##' @importFrom stats setNames
##' 
##' @noRd
.extract_mgf_spectrum <- function(mgf) {
    ## grep description
    desc.idx <- grep("=", mgf)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]

    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+"))
    mode(ms) <- "double"

    if (!length(ms))
        ms <- matrix(numeric(), ncol = 2L)

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)), substring(desc, 1L, r - 1L))
    title <- unname(desc["TITLE"])

    desc[c("PEPMASSMZ", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+")[[1L]][1:2]

    ## select only values of interest and convert to numeric
    desc["CHARGE"] <- sub("[+-]", "", desc["CHARGE"])
    voi <- c("RTINSECONDS", "CHARGE", "SCANS", "PEPMASSMZ", "PEPMASSINT")
    desc <- setNames(as.numeric(desc[voi]), voi)
    desc[is.na(desc[voi])] <- 0L
    list(rtime = unname(desc["RTINSECONDS"]),
         scanIndex = unname(as.integer(desc["SCANS"])),
         precursorMz = unname(desc["PEPMASSMZ"]),
         precursorIntensity = unname(desc["PEPMASSINT"]),
         precursorCharge = unname(as.integer(desc["CHARGE"])),
         mz = ms[, 1L],
         intensity = ms[, 2L],
         title = title)
}
