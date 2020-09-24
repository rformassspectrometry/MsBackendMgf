# MsBackendMgf 0.2

## Changes in 0.2.1

- Adapt to `Spectra::spectraData`/`Spectra::peaksData`.

## Changes in 0.2.0

- Add `export` method to support exporting of data from a `Spectra` object in
  mgf file format.
- Import all fields from an MGF file.
- Rename MGF field names to spectra variables based on `spectraVariableMapping`.
- Support for MGF files with different field names.
