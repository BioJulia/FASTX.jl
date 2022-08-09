```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# FASTA index (FAI files)
FASTX.jl supports FASTA index (FAI) files.
When a FASTA file is indexed with a FAI file, one can seek records by their name, or extract parts of records easily.

See the FAI specifcation here: http://www.htslib.org/doc/faidx.html

### Making an `Index`
A FASTA index (of type `Index`) can be constructed from an `IO` object representing a FAI file:

```jldoctest
julia> io = IOBuffer("seqname\t9\t2\t6\t8");

julia> Index(io) isa Index
true
```

Or from a path representing a FAI file:
```julia
julia> Index("/path/to/file.fai")
```

Alternatively, a FASTA file can be indexed to produce an `Index` using `faidx`.

```jldoctest
julia> faidx(IOBuffer(">abc\nTAGA\nTA"))
Index:
  abc	6	5	4	5
```

Alternatively, a FASTA file can be indexed, and the index immediately written to a FAI file,
by passing an `AbstractString` to `faidx`:

```julia
julia> ispath("/path/to/fasta.fna.fai")
false

julia> faidx("/path/to/fasta.fna");

julia> ispath("/path/to/fasta.fna.fai")
true
```

Note that the restrictions on FASTA files for indexing are stricter than Julia's FASTA parser,
so not all FASTA files that can be read can be indexed:

```jldoctest
julia> str = ">\0\n\0";

julia> first(FASTAReader(IOBuffer(str))) isa FASTARecord
true

julia> Index(IOBuffer(str))
ERROR
[...]
```

### Attaching an `Index` to a `Reader`
When opening a `FASTA.Reader`, you can attach an `Index` by passing the `index` keyword.
You can either pass an `Index` directly, or else an `IO`, in which case an `Index` will be parsed from the `IO`,
or an `AbstractString` that will be interpreted as a path to a FAI file:

```jldoctest
julia> str = ">abc\nTAG\nTA";

julia> idx = faidx(IOBuffer(str));

julia> rdr = FASTAReader(IOBuffer(str), index=idx);
```

### Seeking using an `Index`
With an `Index` attached to a `Reader`, you can do the following operation in O(1) time.
In these examples, we will use the following FASTA file:

```
>seq1 sequence
TAGAAAGCAA
TTAAAC
>seq2 sequence
AACGG
UUGC
```

```@meta
DocTestSetup = quote
using FASTX

data = """>seq1 sequence
TAGAAAGCAA
TTAAAC
>seq2 sequence
AACGG
UUGC
"""

reader = FASTA.Reader(IOBuffer(data), index=faidx(IOBuffer(data)))

end
```

* Seek to a Record using its identifier:
```jldoctest
julia> seekrecord(reader, "seq2");

julia> record = first(reader); sequence(record)
"AACGGUUGC"
```

* Directly extract a record using its identifier
```jldoctest
julia> record = reader["seq1"];

julia> description(record)
"seq1 sequence"
```

* Extract a sequence directly without loading the whole record into memory.
  This is useful for huge sequences like chromosomes
```jldoctest
julia> extract(reader, "seq1", 3:5)
"GAA"
```

```@meta
DocTestSetup = nothing
```

FASTX.jl does not yet support indexing FASTQ files.