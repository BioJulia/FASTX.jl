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
"FASTX" is a shorthand for the two related formats FASTA and FASTQ,
which are handled by the two modules `FASTX.FASTA` and `FASTX.FASTQ`, respectively.

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
```jldoctest
julia> validate_fasta(IOBuffer(">ABC\nDEF")) === nothing
true

julia> validate_fastq(IOBuffer("@ABC\nTAG\n+\nDDD")) === nothing
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
```julia
FASTQWriter(open(path, "w")) do writer
    for record in records
        write(writer, record)
    end
end
```

See more details in the sections in the sidebar.

## Contributing
We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.
