# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `Base.copyto!` methods for copying record data to LongSequences.
- `FASTA.seqlen` & `FASTQ.seqlen` for getting the length of a sequence in a record.

### Changed
- Use BioSequence.jl v2.0 or higher.
- Use TranscodingStreams v0.9.5.

## [1.0.0] - 2019-06-30
### Added
- FASTA submodule.
- FASTQ submodule.
- User manual.
- API reference.

[Unreleased]: https://github.com/BioJulia/FASTX.jl/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/BioJulia/FASTX.jl/tree/v1.0.0
