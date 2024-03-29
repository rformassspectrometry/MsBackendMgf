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

    res <- MsBackendMgf:::.extract_mgf_spectrum(mgf[begin[1]:end[1]])
    expect_true(is.data.frame(res))
    expect_equal(names(res), c("TITLE", "PEPMASS", "CHARGE",
                               "RTINSECONDS", "SCANS", "PEPMASSINT",
                               "mz", "intensity"))
    expect_true(is.na(res$PEPMASSINT))

    res <- MsBackendMgf:::.extract_mgf_spectrum(mgf[begin[2]:end[2]])
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
