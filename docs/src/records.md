```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# Records
FASTX files are considered a sequence of `Record`s, `FASTA.Record` for FASTA files and `FASTQ.Record` for FASTQ.
For convenience, `FASTARecord` and `FASTQRecord` are aliases of `FASTA.Record` and `FASTQ.Record`.

A `Record` object represent the text of the FASTX record as it is, e.g the following FASTA record:
```
>some header here
TAGATGAA
AA
```
Is stored in a `FASTA.Record` object roughly as its constituent bytes, plus some metadata.
There is no notion in the record object of being a DNA or RNA sequence - it's simply an array of bytes.

Records can be constructed from raw parts (i.e. description and sequence and, for FASTQ, quality), where
* `description::AbstractString`
* `sequence::Union{AbstractString, BioSequence}`
* `quality::Union{AbstractString, Vector{<:Number}}`

Alternatively, they can be parsed directly from a string or an `AbstractVector{UInt8}`.

```jldoctest
julia> record = parse(FASTARecord, ">abc\nAGCC\nCCGA");

julia> record2 = FASTARecord("abc", "AGCCCCGA");

julia> record == record2
true
```

Records can be queried for their information, namely identifier, description and sequence (and quality, for FASTQ).
By default, this returns an `AbstractString` view into the `Record`'s data:
```jldoctest
julia> record = parse(FASTARecord, ">ident desc\nUGU\nGA");

julia> (identifier(record), description(record), sequence(record))
("ident", "ident desc", "UGUGA")
```

However, you can ask for getting the sequences as a `String` or any subtype of `BioSequence`:
```jldoctest
julia> record = parse(FASTARecord, ">abc\nUGC\nCCA");

julia> using BioSequences # LongRNA defined in BioSequences.jl

julia> sequence(LongRNA{2}, record)
6nt RNA Sequence:
UGCCCA

julia> sequence(String, record)
"UGCCCA"
```

The number of bytes in the sequence of a `Record` can be queried using `seqsize`:
```jldoctest
julia> record = parse(FASTARecord, ">abc\nUGC\nCCA");

julia> seqsize(record)
6
```

## Byte-oriented sequence positions

FASTX stores and addresses record sequences as raw bytes. It assumes sequence
symbols are ASCII; it does not interpret sequence text as UTF-8. Therefore, the
optional `part` passed to `sequence`, the source offset passed to `copyto!` when
copying from a record to a `BioSequence`, and the range passed to indexed FASTA
`extract` are all one-based byte positions. A byte range can split a multi-byte
UTF-8 symbol and yield invalid text.

If your data deliberately uses multi-byte text symbols, first materialize the
record sequence as a `String`, then build or index a character-oriented
sequence from that string:

```julia
text = sequence(String, record)
characters = collect(text)
characters[2:3]
```

`FASTA.Writer(width=...)` likewise wraps after `width` bytes, rather than after
`width` Unicode characters.

### Reference:
```@docs
identifier
description
sequence
seqsize
```
