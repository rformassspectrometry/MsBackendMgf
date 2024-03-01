# MsBackendMgf 1.11

## Changes in 1.11.2

- Import generic methods from `ProtGenerics`. Requires `ProtGenerics` version
  1.35.3.

## Changes in 1.11.1

- Small runtime improvements in MGF importer.

## Changes in 1.11.0

- Bioconductor 3.19 developmental branch.


# MsBackendMgf 1.7

## Changes in 1.7.5

- Support import of MS levels from MGF files.

## Changes in 1.7.4

- Fix bug in import when m/z values are not sorted.

## Changes in 1.7.3

- Enforce m/z values to be increasingly ordered while importing.

## Changes in 1.7.2

- Run the full unit test suite from the `Spectra` package to check validity of
  the `MsBackendMgf`.


## Changes in 1.7.1

- Use `compress = FALSE` for `NumericList`.

## Changes in 1.7.0

- Bioconductor 3.17 developmental branch.


# MsBackendMgf 1.5

## Changes in 1.5.1

- Fix `export` method to fail if one or more columns contain either S4 classes
  or `list`-like structures.

# MsBackendMgf 1.3

## Changes in 1.3.3

- Import `coreSpectraVariables` from `Spectra`.

## Changes in 1.3.2

- Adapt the `spectraVariableMapping` to `Spectra` version >= 1.5.8.

## Changes in 1.3.1

- Run additional unit test suits defined in `Spectra`.

# MsBackendMgf 1.1

## Changes in 1.1.3

- Fix issue with MGF files lacking peak data.

## Changes in 1.1.2

- Export precursor charge in the expected format (issue
  [#16](https://github.com/rformassspectrometry/MsBackendMgf/issues/16)).

## Changes in 1.1.1

- Add an example to the vignette describing how to export only selected spectra
  variables to the MGF file.

# MsBackendMgf 0.99

## Changes in 0.99.3

- Export `readMgf` function.

## Changes in 0.99.2

- Update installation description in the vignette.

## Changes in 0.99.1

- Directly call internal function from `MsBackendMgf` to avoid parallel
  processing error (function not found) on Windows.


# MsBackendMgf 0.2

## Changes in 0.2.2

- Add examples how to use the `MsBackendMgf` backend to export data to mgf
  files (issue
  [#4](https://github.com/rformassspectrometry/MsBackendMgf/issues/4).


## Changes in 0.2.0

- Add `export` method to support exporting of data from a `Spectra` object in
  mgf file format.
- Import all fields from an MGF file.
- Rename MGF field names to spectra variables based on `spectraVariableMapping`.
- Support for MGF files with different field names.
