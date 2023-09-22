```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# FASTX formatted files

### Readers and writers
A `Reader` and a `Writer` are structs that wrap an IO, and allows efficient reading/writing of FASTX `Record`s.
For FASTA, use `FASTAReader` and `FASTAWriter`, and for FASTQ - well I'm sure you've guessed it.

Readers and writers take control over the underlying IO, and manipulating the IO underneath a Reader/Writer, e.g. by flushing or closing it, cause them to behave in an undefined manner.
Instead, you must interact directly with the reader and writer object.

Closing readers/writers closes the underlying IO.
Because they carry their own buffers, it's important to remember to close writers in particular, else the results may not be fully written to the file.

Readers are iterables of `Record`:

```jldoctest
julia> reader = FASTAReader(IOBuffer(">A\nTAG\n>B\nAGA"));

julia> record = first(reader); typeof(record) == FASTA.Record
true

julia> sequence(record)
"TAG"

julia> # NB! Readers are mutable iterators as can be seen here:

julia> sequence(first(reader))
"AGA"

julia> iterate(reader) === nothing # now empty
true

julia> close(reader)
```

They are normally more than fast enough as they are.
To squeeze extra performance out, you can pass the keyword `copy=false`.
This will cause the reader to return the _same_ record over and over, and mutate it into place.

```jldoctest
julia> reader = FASTAReader(IOBuffer(">A\nTAG\n>B\nAGA"); copy=false);

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

When using readers and writers, be careful that they carry their own buffer,
meaning that the underlying IO may not be updated immediately after reading/writing:
```jldoctest
julia> io = IOBuffer();

julia> writer = FASTAWriter(io);

julia> write(writer, parse(FASTARecord, ">ABC\nDEF"));

julia> take!(io) # NB: Empty!
UInt8[]
```

To use it correctly, either call `flush`, or close the writer first (which also closes the underlying stream).
It is recommended to use readers and writers to `do` syntax in the form:
```jldoctest
julia> FASTAWriter(open(tempname(), "w")) do writer
           write(writer, FASTARecord("identifier", "seq"))
       end
16
```

Which will work for most underlying IO types, and will close the writer when the function returns (hence also closing the underlying IO).

Alternatively, the following syntax may be used:
```jldoctest
julia> open(FASTAWriter, tempname()) do writer
           write(writer, FASTARecord("identifier", "seq"))
       end
16
```

However, this latter syntax does not easily extend to different types of IO, such as gzip compressed streams.

### Validate files
The functions `validate_fasta` and `validate_fastq` can be used to check if an `IO`
contains data that can be read as FASTX.
They return `nothing` if the IO is correctly formatted, and another value if not.

They are significantly faster than parsing the whole file into records,
and are memory efficient.
Be aware that the validators mutate the IO by reading it, so make sure to reset the IO before using it to parse FASTX files.

```jldoctest
julia> io = IOBuffer(">header\r\nAGG\nKK");

julia> validate_fasta(io) === nothing
true

julia> read(io) # NB: IO is now exhausted
UInt8[]

julia> validate_fastq(IOBuffer("@header\nTTT\n+\njkm")) === nothing
true
```

However, this latter syntax does not easily extend to different types of IO, such as gzip compressed streams.

### Validate files
The functions `validate_fasta` and `validate_fastq` can be used to check if an `IO`
contains data that can be read as FASTX.
They return `nothing` if the IO is correctly formatted, and another value if not.

They are significantly faster than parsing the whole file into records,
and are memory efficient.
Be aware that the validators mutate the IO by reading it, so make sure to reset the IO before using it to parse FASTX files.

```jldoctest
julia> io = IOBuffer(">header\r\nAGG\nKK");

julia> validate_fasta(io) === nothing
true

julia> read(io) # NB: IO is now exhausted
UInt8[]

julia> validate_fastq(IOBuffer("@header\nTTT\n+\njkm")) === nothing
true
```
