test_that("readMgf works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")

    expect_error(readMgf(fls), "Please provide a single mgf file.")

    res1 <- readMgf(fls[1])
    cns <- c("rtime", "acquisitionNum", "precursorMz", "precursorIntensity",
             "precursorCharge", "mz", "intensity", "dataOrigin",
             "msLevel", "TITLE")
    expect_identical(sort(names(res1)), sort(cns))
    expect_true(is(res1$mz, "NumericList"))
    expect_true(is(res1$intensity, "NumericList"))
    expect_equal(length(res1$intensity[[1]]), length(res1$mz[[1]]))
    expect_equal(length(res1$intensity[[2]]), length(res1$mz[[2]]))
    expect_equal(length(res1$intensity[[3]]), length(res1$mz[[3]]))
    expect_identical(res1$TITLE,
                     c("File193 Spectrum1719 scans: 2162",
                       "File193 Spectrum1944 scans: 2406",
                       "File193 Spectrum1968 scans: 2432"))

    res2 <- readMgf(fls[2])
    cns <- c("rtime", "precursorMz", "precursorIntensity",
             "precursorCharge", "mz", "intensity", "dataOrigin",
             "msLevel", "TITLE")
    expect_identical(sort(names(res2)), sort(cns))
    expect_true(is(res2$mz, "NumericList"))
    expect_true(is(res2$intensity, "NumericList"))
    expect_equal(length(res2$intensity[[1]]), length(res2$mz[[1]]))
    expect_equal(length(res2$intensity[[2]]), length(res2$mz[[2]]))
    expect_equal(length(res2$intensity[[3]]), length(res2$mz[[3]]))
    expect_true(is.na(res2$rtime[4]))

    res3 <- readMgf(fls[3])
    expect_true(length(res3$mz[[1L]]) == 0)
    expect_true(length(res3$intensity[[1L]]) == 0)
})

test_that("readMgfSplit works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")

    expect_error(readMgfSplit(fls), "Please provide a single mgf file.")
    expect_error(readMgfSplit(fls[1], n = -3), "has to be an integer")
    expect_error(readMgfSplit(fls[1], n = 3), "Consider increasing the value")

    res1 <- readMgfSplit(fls[1])
    cns <- c("rtime", "acquisitionNum", "precursorMz", "precursorIntensity",
             "precursorCharge", "mz", "intensity", "dataOrigin",
             "msLevel", "TITLE")
    expect_identical(sort(names(res1)), sort(cns))
    expect_true(is(res1$mz, "NumericList"))
    expect_true(is(res1$intensity, "NumericList"))
    expect_equal(length(res1$intensity[[1]]), length(res1$mz[[1]]))
    expect_equal(length(res1$intensity[[2]]), length(res1$mz[[2]]))
    expect_equal(length(res1$intensity[[3]]), length(res1$mz[[3]]))
    expect_identical(res1$TITLE,
                     c("File193 Spectrum1719 scans: 2162",
                       "File193 Spectrum1944 scans: 2406",
                       "File193 Spectrum1968 scans: 2432"))
    expect_equal(res1, readMgf(fls[1]))

    res2 <- readMgfSplit(fls[2])
    cns <- c("rtime", "precursorMz", "precursorIntensity",
             "precursorCharge", "mz", "intensity", "dataOrigin",
             "msLevel", "TITLE")
    expect_identical(sort(names(res2)), sort(cns))
    expect_true(is(res2$mz, "NumericList"))
    expect_true(is(res2$intensity, "NumericList"))
    expect_equal(length(res2$intensity[[1]]), length(res2$mz[[1]]))
    expect_equal(length(res2$intensity[[2]]), length(res2$mz[[2]]))
    expect_equal(length(res2$intensity[[3]]), length(res2$mz[[3]]))
    expect_true(is.na(res2$rtime[4]))
    expect_equal(res2, readMgf(fls[2]))

    res3 <- readMgfSplit(fls[3])
    expect_true(length(res3$mz[[1L]]) == 0)
    expect_true(length(res3$intensity[[1L]]) == 0)
    expect_equal(res3, readMgf(fls[3]))
})

test_that(".extract_mgf_spectrum works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")
    mgf <- scan(fls[1], what = "", sep = "\n", quote = "",
                allowEscapes = FALSE, quiet = TRUE)
    cmts <- grep("^[#;!/]", mgf)
    if (length(cmts))
        mgf <- mgf[-cmts]

    begin <- grep("BEGIN IONS", mgf) + 1L
    end <- grep("END IONS", mgf) - 1L
    n <- length(begin)

    res <- .extract_mgf_spectrum(mgf[begin[1]:end[1]])
    expect_true(is.data.frame(res))
    expect_equal(names(res), c("TITLE", "PEPMASS", "CHARGE",
                               "RTINSECONDS", "SCANS", "PEPMASSINT",
                               "mz", "intensity"))
    expect_true(is.na(res$PEPMASSINT))

    with_mocked_bindings(
        ".is_unsorted" = function(...) TRUE,
        code = expect_equal(res, .extract_mgf_spectrum(mgf[begin[1]:end[1]]))
    )

    res <- .extract_mgf_spectrum(mgf[begin[2]:end[2]])
    expect_true(is.data.frame(res))
    expect_equal(names(res), c("TITLE", "PEPMASS", "CHARGE",
                               "RTINSECONDS", "SCANS", "PEPMASSINT",
                               "mz", "intensity"))
    expect_true(!is.na(res$PEPMASSINT))
})

test_that(".format_charge works", {
    res <- .format_charge(c("4+", "3-", NA, "45"))
    expect_equal(res, c("4", "-3", NA_character_, "45"))
})

test_that(".is_unsorted works", {
    x <- cbind(1:3, 1:3)
    expect_false(.is_unsorted(x))

    x <- cbind(c(1, 43, 2), 1:3)
    expect_true(.is_unsorted(x))

    expect_false(.is_unsorted(x[integer(), , drop = FALSE]))
})

test_that(".extract_mgf_spectrum_with_annotations works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")
    a <- grep("spectra.mgf", fls, value = TRUE)
    A <- readLines(a)
    mgf <- A[grep("BEGIN IONS", A)[1L]:grep("END IONS", A)[1L]]
    res <- .extract_mgf_spectrum_with_annotations(mgf)
    expect_true(is.data.frame(res))
    expect_true(length(res$mz) > 0)
    expect_true(length(res$intensity) > 0)
    expect_true(is.data.frame(res$ann_mat[[1L]]))
    expect_true(ncol(res$ann_mat[[1L]]) == 0L)
    expect_equal(nrow(res$ann_mat[[1L]]), length(res$mz[[1L]]))

    b <- grep("fiora", fls, value = TRUE)
    B <- readLines(b)
    mgf <- B[grep("BEGIN IONS", B)[1L]:grep("END IONS", B)[1L]]
    res <- .extract_mgf_spectrum_with_annotations(mgf)
    expect_true(is.data.frame(res))
    expect_true(length(res$mz) > 0)
    expect_true(length(res$intensity) > 0)
    expect_true(length(res$ann_mat) > 0)
    expect_true(is.data.frame(res$ann_mat[[1L]]))
    expect_equal(length(res$intensity[[1L]]), length(res$mz[[1L]]))
    expect_equal(length(res$intensity[[1L]]), nrow(res$ann_mat[[1L]]))

    mgf <- c("BEGIN IONS",
             "TEST=12.3",
             "1 431.1 a b",
             "2 343.1 d",
             "3 313.1",
             "4 542.1 e",
             "END IONS")
    res <- .extract_mgf_spectrum_with_annotations(mgf)
    expect_true(is.data.frame(res))
    expect_true(length(res$mz) > 0)
    expect_true(length(res$intensity) > 0)
    expect_true(length(res$ann_mat) > 0)
    expect_true(is.data.frame(res$ann_mat[[1L]]))
    expect_equal(length(res$intensity[[1L]]), length(res$mz[[1L]]))
    expect_equal(length(res$intensity[[1L]]), nrow(res$ann_mat[[1L]]))
    expect_equal(ncol(res$ann_mat[[1L]]), 2)
})

test_that("readMgf works with annotated = TRUE", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")

    b <- grep("fiora", fls, value = TRUE)
    l <- readLines(b)
    l2 <- c(l, l)
    l2[13] <- "888.575459486 0.002360690850764513"
    tf <- tempfile()
    writeLines(l2, tf)
    res <- readMgf(tf, annotated = TRUE)
    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) == 2)
    expect_equal(res$mz[[1L]], res$mz[[2L]])
    expect_true(is.na(res$V1[[1L]][6]))
    expect_equal(res$V1[[1L]][-6L], res$V1[[2L]][-6L])

    ## Manually defining
    mgf <- c("BEGIN IONS",
             "TITLE=a",
             "1 431.1 a b",
             "2 343.1 d",
             "3 313.1",
             "4 542.1 e",
             "END IONS",
             "BEGIN IONS",
             "TITLE=b",
             "1 2 d",
             "2 2",
             "END IONS",
             "BEGIN IONS",
             "TITLE=c",
             "1 3",
             "2 3",
             "3 3",
             "END IONS"
             )
    writeLines(mgf, tf)
    res <- readMgf(tf, annotated = TRUE)
    expect_equal(res$TITLE, c("a", "b", "c"))
    expect_equal(res$V1, list(c("a", "d", NA, "e"),
                              c("d", NA),
                              c(NA_character_, NA_character_, NA_character_)))
    expect_equal(res$V2, list(c("b", NA, NA, NA),
                              c(NA_character_, NA_character_),
                              c(NA_character_, NA_character_, NA_character_)))
    file.remove(tf)
})
