##' @title Reading MGF files
##'
##' @description
##'
##' The `readMgf` function imports the data from a file in MGF format reading
##' all specified fields and returning the data as a [DataFrame()].
##'
##' @param f `character(1)` with the path to an mgf file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param mapping named `character` vector to rename mgf fields to spectra
##'     variables.
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @return
##'
##' A `DataFrame` with each row containing the data from one spectrum
##' in the MGF file. m/z and intensity values are available in columns `"mz"`
##' and `"intensity"` in a list representation.
##'
##' @export
##'
##' @importFrom Spectra coreSpectraVariables
##'
##' @importFrom S4Vectors DataFrame
##'
##' @importFrom IRanges NumericList
##'
##' @importFrom MsCoreUtils rbindFill
##'
##' @importFrom methods as
##'
##' @author Laurent Gatto, Johannes Rainer
##'
##' @examples
##'
##' fls <- dir(system.file("extdata", package = "MsBackendMgf"),
##'     full.names = TRUE, pattern = "mgf$")[1L]
##'
##' readMgf(fls)
readMgf <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMgf()), ...) {
    requireNamespace("MsBackendMgf", quietly = TRUE)
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
    cmts <- grep("^[#;!/]", mgf, perl = TRUE)
    if (length(cmts))
        mgf <- mgf[-cmts]

    begin <- grep("BEGIN IONS", mgf, fixed = TRUE) + 1L
    end <- grep("END IONS", mgf, fixed = TRUE) - 1L
    n <- length(begin)
    sp <- vector("list", length = n)

    for (i in seq(along = sp))
        sp[[i]] <- .extract_mgf_spectrum(mgf[begin[i]:end[i]],
                                         mapping = mapping)

    res <- DataFrame(rbindFill(sp))

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if(!"msLevel" %in% colnames(res)) {
      res$msLevel <- as.integer(msLevel)
    }

    res
}

readMgf2 <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMgf()),
                    nlines = 100000, BPPARAM = SerialParam(), ...) {
    requireNamespace("MsBackendMgf", quietly = TRUE)
    if (length(f) != 1L)
        stop("Please provide a single mgf file.")
    if (nlines < 1)
        stop("'nlines' has to be an integer > 1.")
    keep_reading <- TRUE
    nskip <- 0L
    sp <- list()
    while (keep_reading) {
        mgf <- scan(file = f, what = "",
                    sep = "\n", quote = "",
                    allowEscapes = FALSE,
                    skip = nskip,
                    nlines = nlines,
                    quiet = TRUE,
                    blank.lines.skip = FALSE)
        message("length(mgf) ", length(mgf), " nlines ",
                nlines, " nskip ", nskip)
        if (length(mgf) < nlines)
            keep_reading <- FALSE

        begin <- grep("BEGIN IONS", mgf) + 1L
        end <- grep("END IONS", mgf) - 1L
        if (!length(end))
            stop("Could not find 'END IONS'. ",
                 "Consider increasing the value for 'nlines'")
        nskip <- nskip + end[length(end)] + 1L
        begin <- begin[seq_along(end)]
        ## This append operation is very costly; need an alternative.
        sp <- c(sp, bpmapply(FUN = function(i, j, mgf, mapping) {
            .extract_mgf_spectrum(mgf[i:j], mapping = mapping)
        }, begin, end, MoreArgs = list(mgf = mgf, mapping = mapping),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM))
    }
    res <- DataFrame(rbindFill(sp))

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if(!"msLevel" %in% colnames(res)) {
      res$msLevel <- as.integer(msLevel)
    }

    res
}

readMgf3 <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMgf()),
                    nlines = 100000, BPPARAM = SerialParam(), ...) {
    requireNamespace("MsBackendMgf", quietly = TRUE)
    if (length(f) != 1L)
        stop("Please provide a single mgf file.")
    if (nlines < 1)
        stop("'nlines' has to be an integer > 1.")
    keep_reading <- TRUE
    nskip <- 0L
    sp <- list()
    while (keep_reading) {
        mgf <- scan(file = f, what = "",
                    sep = "\n", quote = "",
                    allowEscapes = FALSE,
                    skip = nskip,
                    nlines = nlines,
                    quiet = TRUE,
                    blank.lines.skip = FALSE)
        message("length(mgf) ", length(mgf), " nlines ",
                nlines, " nskip ", nskip)
        if (length(mgf) < nlines)
            keep_reading <- FALSE

        begin <- grep("BEGIN IONS", mgf) + 1L
        end <- grep("END IONS", mgf) - 1L
        if (!length(end))
            stop("Could not find 'END IONS'. ",
                 "Consider increasing the value for 'nlines'")
        nskip <- nskip + end[length(end)] + 1L
        begin <- begin[seq_along(end)]
        ## This append operation is very costly; need an alternative.
        sp[[(length(sp) + 1L)]] <- bpmapply(FUN = function(i, j, mgf, mapping) {
            .extract_mgf_spectrum(mgf[i:j], mapping = mapping)
        }, begin, end, MoreArgs = list(mgf = mgf, mapping = mapping),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    }
    res <- DataFrame(rbindFill(unlist(sp, recursive = FALSE, use.names = FALSE)))

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if(!"msLevel" %in% colnames(res)) {
      res$msLevel <- as.integer(msLevel)
    }

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
.extract_mgf_spectrum <- function(mgf, mapping) {
    ## grep description
    desc.idx <- grep("=", mgf, fixed = TRUE)
    desc <- mgf[desc.idx]
    spec <- mgf[-desc.idx]

    ms <- do.call(rbind, strsplit(spec, "[[:space:]]+", perl = TRUE))
    mode(ms) <- "double"

    if (!length(ms) || length(ms) == 1L)
        ms <- matrix(numeric(), ncol = 2L)

    if(nrow(ms) > 1) {
      if (is.unsorted(ms[, 1L])) {
        ms <- ms[order(ms[, 1L]), , drop = FALSE]
      }
    }

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)),
                     substring(desc, 1L, r - 1L))
    title <- unname(desc["TITLE"])

    desc[c("PEPMASS", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+", perl = TRUE)[[1L]][c(1L, 2L)]

    ## Use all fields in the MGF renaming the ones specified by mapping.
    if ("CHARGE" %in% names(desc))
        desc["CHARGE"] <- .format_charge(desc["CHARGE"])
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

#' Format MGF charge string into an integer compatible format.
#'
#' @param x `character`
#' @return `character`, charge without +/- at the end but - as prefix if needed
#' @noRd
.format_charge <- function(x) {
    res <- sub("[+-]", "", x, perl = TRUE)
    negs <- which(endsWith(x, "-"))
    if (length(negs))
        res[negs] <- paste0("-", res[negs])
    res
}

#' @description
#'
#' Function to export a `Spectra` object in MGF format to `con`.
#'
#' @param x `Spectra`
#'
#' @param con output file.
#'
#' @param mapping named `character` vector that maps from `spectraVariables`
#'    (i.e. `names(mapping)`) to the variable name that should be used in the
#'    MGF file.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraVariables spectraNames peaksData spectraData
#'
#' @noRd
#'
#' @examples
#'
#' spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
#' spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
#' spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
#'
#' sps <- Spectra(spd)
#'
#' .export_mgf(sps)
.export_mgf <- function(x, con = stdout(),
                        mapping = spectraVariableMapping(MsBackendMgf()),
                        exportTitle = TRUE) {
    spv <- spectraVariables(x)
    spd <- spectraData(x, spv[!(spv %in% c("dataOrigin", "dataStorage"))])
    col_not_ok <- !vapply(spd, function(z) is.vector(z) & !is.list(z),
                          logical(1))
    if (any(col_not_ok))
        stop("Column(s) ", paste(colnames(spd)[col_not_ok], collapse = ", "),
             " contain multiple elements per row. Please either drop this ",
             "column or reduce its elements to a single value per row.")
    idx <- match(colnames(spd), names(mapping))
    colnames(spd)[!is.na(idx)] <- mapping[idx[!is.na(idx)]]
    if (any(colnames(spd) == "CHARGE")) {
        sign_char <- ifelse(spd$CHARGE > 0, "+", "-")
        nas <- is.na(spd$CHARGE)
        spd$CHARGE <- paste0(abs(spd$CHARGE), sign_char)
        spd$CHARGE[nas] <- ""
    }
    if (!exportTitle)
        spd$TITLE <- NULL
    l <- nrow(spd)
    tmp <- lapply(colnames(spd), function(z) {
        paste0(z, "=", spd[, z], "\n")
    })
    if (exportTitle && !any(colnames(spd) == "TITLE")) {
        if (!is.null(spectraNames(x)))
            title <- paste0("TITLE=", spectraNames(x), "\n")
        else
            title <- paste0("TITLE=msLevel ", spd$msLevel, "; retentionTime ",
                            spd$rtime, "; scanNum ", spd$acquisitionNum, "\n")
        tmp <- c(list(title), tmp)
    }
    pks <- vapply(peaksData(x), function(z)
        paste0(paste0(z[, 1], " ", z[, 2], "\n"), collapse = ""),
        character(1))
    tmp <- do.call(cbind, c(list(rep_len("BEGIN IONS\n", l)),
                            tmp, list(pks), list(rep_len("END IONS\n", l))))
    tmp[grep("=NA\n", tmp)] <- ""
    writeLines(apply(tmp, 1, paste0, collapse = ""), con = con)
}
