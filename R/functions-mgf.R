##' @title Reading MGF files
##'
##' @description
##'
##' The `readMgf()` function imports the data from a file in MGF format reading
##' all specified fields and returning the data as a [S4Vectors::DataFrame()].
##'
##' For very large MGF files the `readMgfSplit()` function might be used
##' instead. In contrast to the `readMgf()` functions, `readMgfSplit()` reads
##' only `nlines` lines from an MGF file at once reducing thus the memory
##' demand (at the cost of a lower performance, compared to `readMgf()`).
##'
##' @param f `character(1)` with the path to an mgf file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param mapping named `character` vector to rename mgf fields to spectra
##'     variables.
##'
##' @param annotated For `readMgf()`: `logical(1)` whether the MGF file
##'     contains additional peak annotations. See examples below or the
##'     documentation for [MsBackendAnnotatedMgf()] for information on the
##'     expected format.
##'
##' @param nlines for `readMgfSplit()`: `integer(1)` with the number of lines
##'     that should be imported and parsed in each iteration.
##'
##' @param BPPARAM parallel processing setup that should be used. Only the
##'     parsing of the imported MGF file is performed in parallel.
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @return
##'
##' A `DataFrame` with each row containing the data from one spectrum
##' in the MGF file. m/z and intensity values are available in columns `"mz"`
##' and `"intensity"` in a list representation. For `readMgf()` with
##' `annotated = TRUE` also all peaks annotation columns (named `"V1", etc)
##' are provided in a list representation, with the lengths of elements
##' matching those of `"mz"` or `"intensity"`.
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
##' @importFrom BiocParallel SerialParam bpmapply
##'
##' @author Laurent Gatto, Johannes Rainer, Sebastian Gibb, Corey Broeckling
##'
##' @examples
##'
##' fls <- dir(system.file("extdata", package = "MsBackendMgf"),
##'     full.names = TRUE, pattern = "mgf$")[1L]
##'
##' readMgf(fls)
##'
##' ## Annotated MGF
##' fl <- system.file("extdata", "xfiora.mgf", package = "MsBackendMgf")
##' res <- readMgf(fl, annotated = TRUE)
##' colnames(res)
##' res$V1
readMgf <- function(f, msLevel = 2L,
                    mapping = spectraVariableMapping(MsBackendMgf()),
                    annotated = FALSE, ...,
                    BPPARAM = SerialParam()) {
    requireNamespace("MsBackendMgf", quietly = TRUE)
    if (length(f) != 1L)
        stop("Please provide a single mgf file.")
    ## Note: using readLines instead has some performance advantages
    ## (few seconds) for very large files
    mgf <- scan(file = f, what = "",
                sep = "\n", quote = "",
                allowEscapes = FALSE,
                quiet = TRUE)

    ## From http://www.matrixscience.com/help/data_file_help.html#GEN
    ## Comment lines beginning with one of the symbols #;!/ can be
    ## included, but only outside of the BEGIN IONS and END IONS
    ## statements that delimit an MS/MS dataset.

    begin <- grep("BEGIN IONS", mgf, fixed = TRUE) + 1L
    end <- grep("END IONS", mgf, fixed = TRUE) - 1L

    if (annotated) fun_extract_mgf <- .extract_mgf_spectrum_with_annotations
    else fun_extract_mgf <- .extract_mgf_spectrum

    res <- bpmapply(begin, end, FUN = function(b, e, mgf)
        fun_extract_mgf(mgf[b:e]), MoreArgs = list(mgf = mgf),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    res <- rbindFill(res)
    if (annotated) {
        p <- as.factor(rep(seq_len(nrow(res)), lengths(res$mz)))
        anns <- rbindFill(res$ann_mat)
        res$ann_mat <- NULL
        pcol <- colnames(anns)
        for (a in pcol)
            res <- do.call(
                "$<-", list(res, name = a, value = unname(split(anns[[a]], p))))
    } else pcol <- character()

    if ("CHARGE" %in% colnames(res))
        res$CHARGE <- .format_charge(res$CHARGE)

    idx <- match(colnames(res), mapping)
    not_na <- !is.na(idx)
    if (any(not_na))
        colnames(res)[not_na] <- names(mapping)[idx][not_na]

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res <- as(res, "DataFrame")
    if (annotated)
        res@metadata <- list(pcol)
    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    if(!"msLevel" %in% colnames(res))
      res$msLevel <- as.integer(msLevel)

    res
}

#' @export
#'
#' @rdname readMgf
readMgfSplit <- function(f, msLevel = 2L,
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
        if (length(mgf) < nlines)
            keep_reading <- FALSE

        begin <- grep("BEGIN IONS", mgf, fixed = TRUE) + 1L
        end <- grep("END IONS", mgf, fixed = TRUE) - 1L
        if (!length(end))
            stop("Could not find 'END IONS'. ",
                 "Consider increasing the value for 'nlines'")
        nskip <- nskip + end[length(end)] + 1L
        begin <- begin[seq_along(end)]
        ## This append operation is very costly; need an alternative.
        sp <- c(sp,
                bplapply(
                    seq_along(begin), function(i, mgf, begin, end) {
                        .extract_mgf_spectrum(mgf[begin[i]:end[i]])},
                    mgf = mgf, begin = begin, end = end,
                    BPPARAM = BPPARAM)
                )
    }
    res <- rbindFill(sp)

    if ("CHARGE" %in% colnames(res))
        res$CHARGE <- .format_charge(res$CHARGE)

    idx <- match(colnames(res), mapping)
    not_na <- !is.na(idx)
    if (any(not_na))
        colnames(res)[not_na] <- names(mapping)[idx][not_na]

    spv <- coreSpectraVariables()
    spv <- spv[!names(spv) %in% c("mz", "intensity")]
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
        if (any(col <- names(spv) == colnames(res)[i]))
            res[[i]] <- as(res[[i]], spv[col][1])
    }

    res <- as(res, "DataFrame")
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
##' @author Laurent Gatto, Johannes Rainer
##'
##' @importFrom stats setNames
##'
##' @noRd
.extract_mgf_spectrum <- function(mgf) {
    ## grep description
    desc.idx <- grep("=", mgf, fixed = TRUE)
    desc <- mgf[desc.idx]

    spec <- strsplit(mgf[-desc.idx], "[[:space:]]+", perl = TRUE)
    if (!length(spec) || length(spec[[1L]]) == 1L)
        ms <- matrix(numeric(), ncol = 2L)
    else
        ms <- matrix(
            as.double(unlist(spec, use.names = FALSE, recursive = FALSE)),
            ncol = length(spec[[1L]]), byrow = TRUE)

    if(.is_unsorted(ms))
        ms <- ms[order(ms[, 1L]), , drop = FALSE]

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)),
                     substring(desc, 1L, r - 1L))

    desc[c("PEPMASS", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+", perl = TRUE)[[1L]][c(1L, 2L)]

    res <- as.data.frame.matrix(matrix(desc, nrow = 1,
                                       dimnames = list(NULL, names(desc))))
    res$mz <- list(ms[, 1L])
    res$intensity <- list(ms[, 2L])
    res
}

#' @title Process MGF files with peak annotations
#'
#' @description
#'
#' Import MS data from MGF files that provide also annotations for individual
#' mass peaks.
#'
#' Expected format of the peak information:
#'
#' - lines with peak information are expected to start with a number (no white
#'   space)
#' - each line is expected to represent one mass peak
#' - the first two elements per line are expected to be the m/z and intensity
#'   values
#' - all following elements represent peak annotation/metadata. For multiple
#'   spectra, they are assumed to be in the same order, i.e., if 3 annotations
#'   are provided for a peak, they are expected to be in the order
#'   annotation_1<white space>annotation_2<white space>annotation_3 for all
#'   spectra. And it is expected that all are provided.
#' - annotations are interpreted as character strings.
#'
#' @author Johannes Rainer and Corey Broeckling
#'
#' @noRd
.extract_mgf_spectrum_with_annotations <- function(mgf) {
    mz <- numeric()
    int <- numeric()
    ann <- as.data.frame(matrix(NA_character_, ncol = 0, nrow = 0))
    ## find peaks
    pks_idx <- grep("^[0-9]", mgf)
    l <- length(pks_idx)
    if (l) {
        pks <- strsplit(mgf[pks_idx], "[[:space:]]+", perl = TRUE)
        ls <- lengths(pks)
        ml <- max(ls)
        if (ml > 2) {
            mli <- 3:ml
            ann <- as.data.frame(
                do.call(rbind, lapply(pks, function(z) z[mli])))
        } else
            ann <- as.data.frame(matrix(NA_character_, ncol = 0, nrow = l))
        mz <- as.numeric(vapply(pks, `[`, NA_character_, 1L, USE.NAMES = FALSE))
        int <-as.numeric(vapply(pks, `[`, NA_character_, 2L, USE.NAMES = FALSE))
        if (is.unsorted(mz)) {
            idx <- order(mz)
            mz <- mz[idx]
            int <- int[idx]
        }
        mgf <- mgf[-pks_idx]
    }
    ## grep description
    desc <- grep("=", mgf, fixed = TRUE, value = TRUE)

    r <- regexpr("=", desc, fixed = TRUE)
    desc <- setNames(substring(desc, r + 1L, nchar(desc)),
                     substring(desc, 1L, r - 1L))

    desc[c("PEPMASS", "PEPMASSINT")] <-
        strsplit(desc["PEPMASS"], "[[:space:]]+", perl = TRUE)[[1L]][c(1L, 2L)]

    res <- as.data.frame.matrix(matrix(desc, nrow = 1,
                                       dimnames = list(NULL, names(desc))))
    res$mz <- list(mz)
    res$intensity <- list(int)
    res$ann_mat <- list(ann)
    res
}

.is_unsorted <- function(x) {
    nrow(x) && is.unsorted(x[, 1L])
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
#' @importMethodsFrom Spectra spectraVariables spectraNames spectraData
#'
#' @importMethodsFrom ProtGenerics peaksData
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
    if (!exportTitle && any(colnames(spd) %in% "TITLE"))
        spd$TITLE <- NULL
    l <- nrow(spd)
    tmp <- lapply(colnames(spd), function(z) {
        paste0(z, "=", spd[, z], "\n")
    })
    if (exportTitle && !any(colnames(spd) %in% "TITLE")) {
        sn <- spectraNames(x)
        if (any(sn != as.character(seq_along(x))))
            title <- paste0("TITLE=", sn, "\n")
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
