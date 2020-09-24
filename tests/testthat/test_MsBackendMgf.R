test_that("backendInitialize,MsBackendMgf works", {
    fls <- dir(system.file("extdata", package = "MsBackendMgf"),
               full.names = TRUE, pattern = "mgf$")
    be <- MsBackendMgf()

    ## Import a single file.
    res1 <- backendInitialize(be, fls[1])
    n1 <- length(res1) ## 3
    expect_identical(length(res1), n1)
    expect_identical(res1$dataStorage, rep("<memory>", n1))
    expect_identical(res1$dataOrigin, rep(normalizePath(fls[1]), n1))
    expect_identical(res1$msLevel, rep(2L, n1))
    expect_identical(lengths(res1$mz), c(14L, 21L, 14L))

    res2 <- backendInitialize(be, fls[2])
    n2 <- length(res2) ## 4

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

    ## TODO: Import with failing file.
    ## TODO: Import with failing file and nonStop = TRUE

    ## errors
    expect_error(backendInitialize(be), "'files' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(backendInitialize(be, "a"), "a not found")
})

test_that("spectraVariableMapping works", {
    res <- spectraVariableMapping()
    expect_true(is.character(res))
    expect_true(length(names(res)) == length(res))

    expect_error(spectraVariableMapping(format = "other"))
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
})
