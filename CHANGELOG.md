# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [2.1.0]
### Additions
* Implement `Base.copy!` for `FASTQRecord` and `FASTARecord`

## [2.0.1]
### Bugfix
* Fix `Base.read!(::FASTQReader, ::FASTQRecord)` (issue #95)

## [2.0.0]
Version 2 is a near-complete rewrite of FASTX.
It brings strives to provide an easier and more consistent API, while also being
faster, more memory efficient, and better tested than v1.

The changes are comprehensive, but code should only need a few minor tweaks to
work with v2. I recommend upgrading your packages using a static analysis tool like JET.jl.

### Breaking changes
#### Records
* `description` has changed meaning: In v1, it meant the part of the header after the '>' symbol
  and up until first whitespace. Now it extends to the whole header line until the ending newline.
  This implies the identifier is a prefix of the description.
* `header` has been removed, and is now replaced by `description`.
* All `Record` objects now have an identifier, a description and a sequence, and all `FASTQRecord`s
  have a quality. These may be empty, but will not throw an error when accessing them.
* As a consequence, all "checker" functions like `hassequence`, `isfilled`, `hasdescription` and
  so on has been removed, since the answer now is trivially "yes" in all cases.
* `identifier`, `description`, `sequence` and `quality` now returns an `AbstractString` by default.
  Although it is an implementation detail, it uses zero-copy string views for performance.
* You can no longer construct a record using e.g. `Record(::String)`. Instead, use `parse(Record, ::String)`.
* `seqlen` is renamed `seqsize` to emphasize that it returns the data size of the sequence,
  not necessarily its length.

#### Readers/writers
* All readers/writers now take any other arguments than the main IO as a keyword for clarity
  and consistency.
* FASTQ.Writers will no longer by default modify `FASTQ.Records`'s second header.
  An optional keyword forces the reader to always write/skip second header if set to `true` or `false`,
  but it defaults to `nothing`, meaning it leaves it intact.
* FASTQ writers now can no longer fill in ambiguous bases in Records transparently,
  or otherwise transform records, when writing.
  If the user wishes to transform records, they must do it my manually calling a function that transforms the records.

#### Other breaking changes
* `FASTQ.Read` has been removed. To subset a read, extract the sequence and quality, and construct
  a new Record object from these.
* `transcribe` has been removed, as it is now trivial to do the same thing.
  It may be added in a future release with new functionality.

### New features
* Function `quality_scores` return the qualities of a FASTQ record as a lazy, validating iterator
  of PHRED quality scores.
* New object: `QualityEncoding` can be used to construct custom PHRED/ASCII quality encodings.
  accessing quality scores uses an existing default object.
* Readers now have a keyword `copy` that defaults to `true`. If set to `false`, iterating over
  a reader will overwrite the same record for performance. Use with care.
  This makes the old `while !eof(reader)`-idiom obsolete in favor of iterating over a reader
  constructed with `copy=false`.
* Users can now use the following syntax to make processing gzipped readers easier:
  ```
  Reader(GzipDecompressorStream(open(path)); kwargs...) do reader
      # stuff
  end
  ```
  this is a change in BioGenerics.jl, but is guaranteed to work in FASTX.jl v2.
* FAI (FASTX index) files can now be written as well as read.
* FASTA files can now be indexed with the new function `faidx`.
* Function `extract` can extract parts of a sequence from an indexed FASTA reader
  without loading the entire sequence into memory.
  You can use this to e.g. extract a small part of a large chromosome. (see #29)
* New functions `validate_fasta` and `validate_fastq` validates if an `IO` is formatted
  validly, faster and more memory-efficiently than loading in the file.

### Other changes
* All practically useful functions and types are now exported directly from FASTX,
  so users don't need to prepend identifiers with `FASTA.` or `FASTQ.`.
* FASTA readers are more liberal in what formats they will accept (#73)  

### Removed
* The method `FASTA.sequence(::FASTA.Record)` has been removed, since the auto-detection of sequence
  type chould not be made reliable enough.

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
[4;1386;2550t]
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
