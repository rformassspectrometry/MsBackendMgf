library(testthat)
library(MsBackendMgf)

test_check("MsBackendMgf")

## Run tests defined in test suites from the Spectra package.
fls <- dir(system.file("extdata", package = "MsBackendMgf"),
           full.names = TRUE, pattern = "mgf$")[1:2]
be <- MsBackendMgf()
be <- backendInitialize(be, fls)

## Run the MsBackend spectra variable test suite
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")

## Run single test file.
test_file(paste0(test_suite, "/test_spectra_variables.R"))
## Run the whole suite.
## test_dir(test_suite)
