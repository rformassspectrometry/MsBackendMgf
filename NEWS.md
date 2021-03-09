# MsBackendMgf 0.99

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
