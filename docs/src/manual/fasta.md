```@meta
CurrentModule = FASTX
```

# FASTA formatted files

FASTA is a text-based file format for representing biological sequences.
A FASTA file stores a list of sequence records with name, description, and
sequence.

The template of a sequence record is:

```
>{name} {description}?
{sequence}
```

Here is an example of a chromosomal sequence:

```
>chrI chromosome 1
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC
CACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG
```

## Readers and Writers
The reader and writer for FASTA formatted files, are found within the
FASTA submodule of FASTX.

- [`FASTA.Reader`](@ref)
- [`FASTA.Writer`](@ref)

They can be created with IOStreams.

```jlcon
using FASTX

r = FASTA.Reader(open("my-seqs.fasta", "r"))
w = FASTA.Writer(open("my-out.fasta", "w"))
```

Alternatively, `Base.open` is overloaded with a method for conveinience:

```jlcon
r = open(FASTA.Reader, "my-seqs.fasta")
w = open(FASTA.Writer, "my-out.fasta")
```

Usually sequence records will be read sequentially from a file by iteration.

```jlcon
reader = open(FASTA.Reader, "my-seqs.fasta")
for record in reader
    ## Do something
    # like showing the identifiers
    @show FASTA.identifier(record)
end
close(reader)
```

You can also overwrite records in a while loop to avoid excessive memory allocation.

```jlcon
reader = open(FASTA.Reader, "my-seqs.fasta")
record = FASTA.Record()
while !eof(reader)
    read!(reader, record)
    ## Do something.
end
```

But if the FASTA file has an auxiliary index file formatted in fai, the reader
supports random access to FASTA records, which would be useful when accessing
specific parts of a huge genome sequence:

```jlcon
reader = open(FASTA.Reader, "sacCer.fa", index = "sacCer.fa.fai")
chrIV = reader["chrIV"]  # directly read sequences called chrIV.
close(reader)
```

Reading in a sequence from a FASTA formatted file will give you a variable of
type [`FASTA.Record`](@ref).

Various getters and setters are available for [`FASTA.Record`](@ref)s:

- [`FASTA.hasidentifier`](@ref)
- [`FASTA.identifier`](@ref)
- [`FASTA.hasdescription`](@ref)
- [`FASTA.description`](@ref)
- [`FASTA.hassequence`](@ref)
- [`FASTA.sequence`](@ref)

To write a `BioSequence` to FASTA file, you first have to create a [`FASTA.Record`](@ref):

```jlcon
using BioSequences
x = dna"aaaaatttttcccccggggg"
rec = FASTA.Record("MySeq", x)
w = open(FASTA.Writer, "my-out.fasta")
write(w, rec)
close(w)
```

As always with julia IO types, remember to close your file readers and writer
after you are finished.

Using `open` with a do-block can help ensure you close a stream after you are
finished.

```jlcon
open(FASTA.Reader, "my-reads.fasta") do reader
    for record in reader
        ## Do something
    end
end
```
