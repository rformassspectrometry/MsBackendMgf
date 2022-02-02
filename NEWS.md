# MsBackendMgf 1.3

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
