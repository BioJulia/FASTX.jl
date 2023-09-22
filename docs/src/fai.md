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
```jldoctest
julia> Index("../test/data/test.fasta.fai");
```

Alternatively, a FASTA file can be indexed to produce an `Index` using `faidx`.

```jldoctest
julia> faidx(IOBuffer(">abc\nTAGA\nTA"))
Index:
  abc	6	5	4	5
```

Alternatively, a FASTA file can be indexed, and the index immediately written to a FAI file,
by passing an `AbstractString` to `faidx`:

```jldoctest
julia> rm("../test/data/test.fasta.fai") # remove existing fai

julia> ispath("../test/data/test.fasta.fai")
false

julia> faidx("../test/data/test.fasta");

julia> ispath("../test/data/test.fasta.fai")
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

### Writing a FAI file
If you have an `Index` object, you can simply `write` it to an IO:
```jldoctest
julia> index = open(i -> Index(i), "../test/data/test.fasta.fai");

julia> filename = tempname();

julia> open(i -> write(i, index), filename, "w");

julia> index2 = open(i -> Index(i), filename);

julia> string(index) == string(index2)
true
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

You can also add a index to an existing reader using the `index!` function:

```@docs
index!
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

### Reference:
```@docs
faidx
seekrecord
extract
Index
```
