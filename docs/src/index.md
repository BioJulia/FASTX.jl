# FASTX
[![Latest Release](https://img.shields.io/github/release/BioJulia/FASTX.jl.svg)](https://github.com/BioJulia/FASTX.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/FASTX.jl/blob/master/LICENSE) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3663087.svg)](https://doi.org/10.5281/zenodo.3663087)
[![Pkg Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Chat](https://img.shields.io/gitter/room/BioJulia/FASTX.svg)](https://gitter.im/BioJulia/FASTX.jl)

Read and write files in FASTA and FASTQ format, the most common biological sequence file format.

## Installation
You can install FASTX from the julia REPL.
Press `]` to enter pkg mode again, and enter the following:

```julia
(v1.8) pkg> add FASTX
```

## Quickstart
* Construct records from raw parts
```jldoctest
julia> record = FASTARecord("some header", dna"TAGAAGA");

julia> (identifier(record), description(record), sequence(record))
("some", "some header", "TAGAAGA")

julia> sequence(LongDNA{2}, record)
7nt DNA Sequence:
TAGAAGA
```

* Validate files
```julia
julia> nothing === validate_fasta(IOBuffer(">ABC\nDEF"))
true
```

* Read FASTX files
```julia
record = FASTAReader(first, IOBuffer(">ABC\nDEF"))

sequence(record) == "DEF" # should be true

# Or with do-syntax
FASTAReader(GzipDecompressorStream(open(path))) do reader
    for record in reader
        # do something with record
    end
end
```

* Write FASTX files
```jldoctest
julia> buffer = IOBuffer()

julia> FASTQWriter(buffer) do writer
        write(writer, parse(FASTQRecord, "@header\nTAG\n+\nJJK"))
        flush(writer)
        String(take!(buffer))
    end
"@header\nTAG\n+\nJJK"

```

See more details below

## FASTX overview
### Records
"FASTX" is a shorthand for the two related formats FASTA and FASTQ,
which are handled by the two modules `FASTX.FASTA` and `FASTX.FASTQ`, respectively.

The FASTA and FASTQ modules have very similar design decisions and abstractions, covered in this section.
For convenience, examples here will use FASTA.
For more details on each of the modules, see their dedicated sections in the sidebar.

FASTX files are considered a sequence of `Record`s.
A `Record` object represent the text of the FASTX record as it is, e.g the following FASTA record:
```
>some header here
TAGATGAA
AA
```
Is stored in a `FASTA.Record` object.
Records can be constructed from raw parts (i.e. description and sequence and, for FASTQ, quality), where
* `description::AbstractString`
* `sequence::Union{AbstractString, BioSequence}`
* `quality::Union{AbstractString, Vector{<:Number}}`

Alternatively, they can be parsed directly from a string or an `AbstractVector{UInt8}`.

```jldoctest
julia> record = parse(FASTARecord, ">abc\nAGCC\nCCGA");

julia> record2 = FASTARecord("abc", dna"AGCCCCGA");

julia> record == record2
true
```

Records can be queried for their information, namely identifier, description and sequence (and quality, for FASTQ).
By default, this returns an `AbstractString` view into the `Record`'s data:
```jldoctest
julia> record = parse(FASTARecord, "ident desc\nUGU\nGA");

julia> (identifier(record), description(record), sequence(record))
("ident", "ident desc", "UGUGA")
```

However, you can ask for getting the sequences as a `String` or any subtype of `BioSequence`:
```jldoctest
julia> record = parse(FASTARecord, ">abc\nUGC\nCCA");

julia> sequence(LongRNA{2}, record)
6nt RNA Sequence:
UGCCCA

julia> sequence(String, record)
"UGCCCA"
```

The length of the sequence (in bytes, _not_ in characters!) can be accessed with `seqlen(record)`:
```jldoctest
julia> record = parse(FASTARecord, ">X\ntaga\nccaa");

julia> seqlen(record)
8
```

### Validate files
The functions `validate_fasta` and `validate_fastq` can be used to check if an IO
contains data that can be read as FASTX.
They are significantly faster than parsing the whole file into records,
and are memory efficient.

Be aware that the validators mutate the IO by reading it, so make sure to reset the IO before using it to parse FASTX files.

### Readers and writers
A `Reader` and a `Writer` are structs that wrap an IO, and allows efficient reading/writing of FASTX `Record`s.
Both these types take control over the underlying IO, and manipulating the IO underneath a Reader/Writer cause them to behave in an undefined manner.

Closing them closes the underlying stream.
Because they carry their own buffers, it's important to remember to close writers in particular.

Readers are iterables of `Record`:

```jldoctest
julia> reader = FASTAReader(IOBuffer(">A\nTAG\n>B\nAGA"))

julia> sequence(first(reader))
"TAG"

julia> # NB! They are mutable iterators as can be seen here:

julia> sequence(first(reader))
"AGA"

julia> iterate(reader) === nothing
true
```

They are normally more than fast enough as they are.
To squeeze extra performance out, you can pass the keyword `copy=false`.
This will cause the reader to return the _same_ record over and over, and mutate it into place.

```jldoctest
julia> reader = FASTAReader(IOBuffer(">A\nTAG\n>B\nAGA"); copy=false)

julia> rec1 = first(reader); sequence(rec1)
"TAG"

julia> rec2 = first(reader); sequence(rec2)
"AGA"

julia> rec1 === rec2
true

julia> sequence(rec1)
"AGA"

julia> close(reader)
```

When using writers, be careful that they carry their own buffer:
```julia
julia> buffer = IOBuffer();

julia> writer = FASTAWriter(buffer);

julia> write(writer, parse(FASTARecord, ">ABC\nDEF"))

julia> take!(buffer) # NB: Empty!
UInt8[]
```

To use it correctly, either call `flush`, or close the writer first (which also closes the underlying stream).
It is recommended to use readers and writers to `do` syntax in the form:
```julia
Writer(open(my_file)) do writer
    for record in my_records
        write(writer, record)
    end
end
```

Which will work for any IO.

## Contributing
We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.
