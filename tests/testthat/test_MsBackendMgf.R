test_that("backendInitialize,MsBackendMgf works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")[1:2]
    be <- MsBackendMgf()

    ## Import a single file.
    res1 <- backendInitialize(be, fls[1])
    n1 <- length(res1) ## 3
    expect_identical(length(res1), n1)
    expect_identical(res1$dataStorage, rep("<memory>", n1))
    expect_identical(res1$dataOrigin, rep(normalizePath(fls[1]), n1))
    expect_identical(res1$msLevel, rep(2L, n1))
    expect_identical(lengths(res1$mz), c(14L, 21L, 14L))

    res1_b <- backendInitialize(be, fls[1], nlines = 30)
    expect_equal(res1, res1_b)

    res2 <- backendInitialize(be, fls[2])
    n2 <- length(res2) ## 4
    expect_equal(n2, 4)

    ## Import multiple files.
    res_all <- backendInitialize(be, fls)
    expect_true(length(res_all) == n1 + n2)
    expect_identical(res_all[1]$mz, res1[1]$mz)
    expect_identical(res_all[n1 + 1]$mz, res2[1]$mz)
    expect_true(all(res_all$msLevel == 2L))
    expect_identical(res_all$dataOrigin,
                     c(rep(normalizePath(fls[1]), n1),
                       rep(normalizePath(fls[2]), n2)))
    expect_true(is.integer(res_all@spectraData$msLevel))

    ## errors
    expect_error(backendInitialize(be), "'files' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(backendInitialize(be, "a"), "a not found")
    expect_error(backendInitialize(be, fls[1], nlines = "a"), "integer")
    expect_error(backendInitialize(be, "a"), "not found")
})

test_that("spectraVariableMapping works", {
    res <- spectraVariableMapping(MsBackendMgf())
    expect_true(is.character(res))
    expect_true(length(names(res)) == length(res))

    expect_error(spectraVariableMapping(MsBackendMgf(), format = "other"))
})

test_that("mixed MS level import works", {

  fls <- dir(system.file("extdata", package = "MsBackendMgf"),
             full.names = TRUE, pattern = "mgf$")[4]

  custom_mapping <- c(rtime = "RTINSECONDS",
                      acquisitionNum = "SCANS",
                      precursorMz = "PEPMASS",
                      precursorIntensity = "PEPMASSINT",
                      precursorCharge = "CHARGE",
                      msLevel = "MSLEVEL")

  res <- Spectra(fls,
                 source = MsBackendMgf(),
                 backend = MsBackendDataFrame(),
                 mapping = custom_mapping)

  expect_identical(length(res), 2L)
  expect_identical(res$msLevel, c(1L, 2L))
})

test_that("export,MsBackendMgf works", {
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))

    sps <- Spectra(spd)

    fl <- tempfile()
    export(MsBackendMgf(), sps, file = fl)
    res <- backendInitialize(MsBackendMgf(), fl)
    expect_equal(rtime(res), rtime(sps))
    expect_equal(peaksData(res), peaksData(sps@backend))

    sps$precursorCharge <- c(2L, NA_integer_, -4L)
    export(MsBackendMgf(), sps, file = fl)
    res <- backendInitialize(MsBackendMgf(), fl)
    expect_equal(precursorCharge(res), precursorCharge(sps))

    export(MsBackendMgf(), sps, file = fl, exportTitle = FALSE)
    res <- readLines(fl)
    expect_true(length(grep("TITLE", res)) == 0)

    sps$TITLE <- c("a", "b", "c")
    export(MsBackendMgf(), sps, file = fl, exportTitle = FALSE)
    res <- readLines(fl)
    expect_true(length(grep("TITLE", res)) == 0)

    spectraNames(sps) <- c("d", "e", "f")
    export(MsBackendMgf(), sps, file = fl, exportTitle = FALSE)
    res <- readLines(fl)
    expect_true(length(grep("TITLE", res)) == 0)

    spectraNames(sps) <- NULL
    export(MsBackendMgf(), sps, file = fl, exportTitle = TRUE)
    res <- readLines(fl)
    expect_true(length(grep("TITLE=msLevel", res)) == 0)

    spectraNames(sps) <- paste0("yes", seq_along(sps))
    export(MsBackendMgf(), sps, file = fl, exportTitle = TRUE)
    res <- readLines(fl)
    expect_true(length(grep("TITLE=yes", res)) == 0)

    expect_error(export(MsBackendMgf(), file = fl), "missing")
    expect_error(export(MsBackendMgf(), x = spd, file = fl), "spectra data to")

    sps$lst <- list(1:3, c(4, 5, 2), c("a", "b"))
    expect_error(export(MsBackendMgf(), sps, file = fl), "multiple elements")

    sps$fail <- IRanges::NumericList(sps$mz, compress = FALSE)
    expect_error(export(MsBackendMgf(), sps, file = fl), "multiple elements")
})

test_that("Spectra works with MsBackendMgf", {
    s <- Spectra(fls, source = MsBackendMgf())
    be <- backendInitialize(MsBackendMgf(), fls)
    expect_equal(be$mz, s$mz)
})
