# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2021-07-13
### Added:
* `header(::Union{FASTA.Record, FASTQ.Record})` returns the full header line.
* `sequence_iter(::Union{FASTA.Record, FASTQ.Record})` returns a no-copy iterator over the sequence. If the record is mutated, this iterator will be in an invalid state.
* `quality_iter(::FASTQ.Record)` - same as above, but for PHRED quality.
* New type `FASTQRead` stores the same data as a FASTQ record, but in a Julia native format instead of a ASCII-encoding byte vector. (PR #35)

### Bugfixes
* Allow trailing newlines after last record of FASTA and FASTQ
* Fix parser FSM ambiguity
* Fix off-by-one error in line counting of FASTQ files
* Various small fixes to the internal parsing regex
* Writers are now parametric and buffered for increased writing speed
* Fixed a bug where Windows-style newlines would break the parser

## Unreleased

## [1.1.0] - 2019-08-07
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

[Unreleased]: https://github.com/BioJulia/FASTX.jl/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/BioJulia/FASTX.jl/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/BioJulia/FASTX.jl/tree/v1.0.0
