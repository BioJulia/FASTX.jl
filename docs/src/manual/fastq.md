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

As always with julia IO types, remember to close your file readers and writer
after you are finished.

Using `open` with a do-block can help ensure you close a stream after you are
finished. `Base.open` is overloaded with a method for this purpose.

```jlcon
r = open(FASTQ.Reader, "my-reads.fastq")
w = open(FASTQ.Writer, "my-out.fastq")
```

Note that [`FASTQ.Reader`](@ref) does not support line-wraps within sequence and quality.
Usually sequence records will be read sequentially from a file by iteration.

```jlcon
open(FASTQ.Reader, "my-reads.fastq") do reader
    for record in reader
        ## Do something
    end
end
```

Gzip compressed files can be streamed to the `Reader`
using the [CodecZlib.jl](https://github.com/JuliaIO/CodecZlib.jl) package.

```jlcon
reader = FASTQ.Reader(GzipDecompressorStream(open("my-reads.fastq.gz")))
for record in reader
    ## do something
end
close(reader)
```

You can also overwrite records in a while loop to avoid excessive memory allocation.

```jlcon
open(FASTQ.Reader, "my-reads.fastq") do reader
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        ## Do something.
    end
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

```jlcon
using BioSequences
x = dna"aaaaatttttcccccggggg"
q = fill(1, length(x))
rec = FASTQ.Record("MySeq", x, q)
open(FASTQ.Writer, "my-out.fastq") do
    write(w, rec)
end
```

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

## FASTQ Reads

The FASTQ Record data structure is very close to the FASTQ file-format, and stores all data in a data vector.
The FASTQ Read data structure is better suited for sequence manipulations:

- Identifier and description are Strings.
- Sequence is stored as a `BioSequence` and not just as ASCII characters.
- The quality is stored as raw PHRED-score (Integer), and there is no offset to worry about after the conversion.

A FASTQ Read record also allows for convenient sub-setting using normal range syntax as for Arrays and Strings:

```jlcon
using FASTX

read = FASTQ.Read(first(FASTQ.Reader(open("my-reads.fastq", "r"))))
# FASTX.FASTQ.FASTQRead{BioSequences.DNAAlphabet{4}}:
#    identifier: SEQ_ID
#   description: 
#      sequence: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
#       quality: [0, 6, 6, 9, 7, 7, 7, 7, 9, 9  …  29, 34, 34, 34, 34, 34, 34, 34, 21, 20]

length(read)
# 60

read[3:12]
# FASTX.FASTQ.FASTQRead{BioSequences.DNAAlphabet{4}}:
#    identifier: SEQ_ID
#   description: 
#      sequence: TTTGGGGTTC
#       quality: [6, 9, 7, 7, 7, 7, 9, 9, 9, 10]

read[3:12].sequence
# 10nt DNA Sequence:
# TTTGGGGTTC

read[1:3].quality
# 3-element Array{UInt8,1}:
#  0x00
#  0x06
#  0x06
```

When defining the FASTQ.Read object, the quality can be given as an integer offset or encoding symbol.

```jlcon
rec = first(FASTQ.Reader(open("my-reads.fastq", "r")));
r1 = FASTQ.FASTQ.Read(rec, 33);
r2 = FASTQ.FASTQ.Read(rec, :illumina18);
r1.quality == r2.quality
# true
```

The platform specific quality encodings are described on [wikipedia](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).
