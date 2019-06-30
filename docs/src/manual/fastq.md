```@meta
CurrentModule = FASTX
```

# FASTQ formatted files

FASTQ is a text-based file format for representing DNA sequences along with
qualities for each base.
A FASTQ file stores a list of sequence records in the following format:

```
@{name} {description}?
{sequence}
+
{qualities}
```

Here is an example of one record from a FASTQ file:

```
@FSRRS4401BE7HA
tcagTTAAGATGGGAT
+
###EEEEEEEEE##E#
```

## Readers and Writers

The reader and writer for FASTQ formatted files, are found within the
FASTQ module of FASTX.

- [`FASTQ.Reader`](@ref)
- [`FASTQ.Writer`](@ref)

They can be created with IOStreams:

```jlcon
using FASTX

r = FASTQ.Reader(open("my-reads.fastq", "r"))
w = FASTQ.Writer(open("my-output.fastq", "w"))
```

Alternatively, `Base.open` is overloaded with a method for conveinience:

```jlcon
r = open(FASTQ.Reader, "my-reads.fastq")
w = open(FASTQ.Writer, "my-out.fastq")
```

Note that [`FASTQ.Reader`](@ref) does not support line-wraps within sequence and quality.
Usually sequence records will be read sequentially from a file by iteration.

```jlcon
reader = open(FASTQ.Reader, "my-reads.fastq")
for record in reader
    ## Do something
end
close(reader)
```

You can also overwrite records in a while loop to avoid excessive memory allocation.

```jlcon
reader = open(FASTQ.Reader, "my-reads.fastq")
record = FASTQ.Record()
while !eof(reader)
    read!(reader, record)
    ## Do something.
end
```

Reading in a record from a FASTQ formatted file will give you a variable of
type [`FASTQ.Record`](@ref).

Various getters and setters are available for [`FASTQ.Record`](@ref)s:

- [`FASTQ.hasidentifier`](@ref)
- [`FASTQ.identifier`](@ref)
- [`FASTQ.hasdescription`](@ref)
- [`FASTQ.description`](@ref)
- [`FASTQ.hassequence`](@ref)
- [`FASTQ.sequence`](@ref)
- [`FASTQ.hasquality`](@ref)
- [`FASTQ.quality`](@ref)

To write a `BioSequence` to FASTQ file, you first have to create a [`FASTQ.Record`](@ref):

As always with julia IO types, remember to close your file readers and writer
after you are finished.

Using `open` with a do-block can help ensure you close a stream after you are
finished.

```jlcon
open(FASTQ.Reader, "my-reads.fastq") do reader
    for record in reader
        ## Do something
    end
end
```

## Quality encodings

FASTQ records have a quality string which have platform dependent encodings.
The FASTQ submodule has encoding and decoding support for the following
quality encodings. These can be used with a [`FASTQ.quality`](@ref) method, to
ensure the correct quality score values are extracted from your FASTQ quality
strings. 

- [`FASTQ.SANGER_QUAL_ENCODING`](@ref)
- [`FASTQ.SOLEXA_QUAL_ENCODING`](@ref)
- [`FASTQ.ILLUMINA13_QUAL_ENCODING`](@ref)
- [`FASTQ.ILLUMINA15_QUAL_ENCODING`](@ref)
- [`FASTQ.ILLUMINA18_QUAL_ENCODING`](@ref)

