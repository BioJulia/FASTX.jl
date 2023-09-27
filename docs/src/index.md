```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX, BioSequences
end
```

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
See more documentation in the sections in the sidebar.

### Read FASTA or FASTQ files
It is preferred to use the `do` syntax to automatically close the file when you're done with it:
```jldoctest
julia> FASTAReader(open("../test/data/test.fasta")) do reader
           for record in reader
               println(identifier(record))
           end
       end
abc
```

Alternatively, you can open and close the reader manually:

```jldoctest
julia> reader = FASTAReader(open("../test/data/test.fasta"));

julia> for record in reader
           println(identifier(record))
       end
abc
       
julia> close(reader)
```

### Write FASTA or FASTQ files
```jldoctest
julia> FASTQWriter(open(tempname(), "w")) do writer
           write(writer, FASTQRecord("abc", "TAG", "ABC"))
       end
15
```

### Read and write Gzip compressed FASTA files
```jldoctest
julia> using CodecZlib

julia> FASTAReader(GzipDecompressorStream(open("../test/data/seqs.fna.gz"))) do reader
           for record in reader
               println(identifier(record))
           end
       end
seqa
seqb

julia> FASTQWriter(GzipCompressorStream(open(tempname(), "w"))) do writer
           write(writer, FASTQRecord("header", "sequence", "quality!"))
       end
28
```

For added convenience, you can also use the reader and writer macros `rdr""` and `wtr""`.
These macros use the file extensions to determine the biological sequence reader or writer type, and any file compresion.
To use these macros with the `do`-syntax, you can use `open` as normal. Hence, the above code block can also be written in the following equivalent way:

```jldoctest
julia> using CodecZlib

julia> open(rdr"../test/data/seqs.fna.gz") do reader
           for record in reader
               println(identifier(record))
           end
       end
seqa
seqb
```

### Construct FASTA or FASTQ records from raw parts
```jldoctest
julia> fasta_record = FASTARecord("some header", dna"TAGAAGA");

julia> fastq_record = FASTQRecord("read1", "TAGA", "ABCD");
```

### Validate that a file (or an arbitrary `IO`) is well-formatted
The `validate_fast*` functions return `nothing` if the IO is well formatted
```jldoctest
julia> validate_fasta(IOBuffer(">ABC\nDEF")) === nothing
true

julia> validate_fastq(IOBuffer("@ABC\nTAG\n+\nDDD")) === nothing
true
```

To check if files are well-formatted:
```jldoctest
julia> open(validate_fasta, "../test/data/test.fasta") === nothing
true

julia> open(validate_fasta, "Project.toml") === nothing
false
```

## Contributing
We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.

