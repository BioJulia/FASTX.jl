```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# FASTQ formatted files
__NB: First read the overview in the sidebar__

FASTQ is a text-based file format for representing DNA sequences along with qualities for each base.
A FASTQ file stores a list of sequence records in the following format:

The template of a sequence record is:

```
@{description}
{sequence}
+{description}?
{qualities}
```

Where the "identifier" is the first part of the description up to the first whitespace
(or the entire description if there is no whitespace)

The description may optionally be present on the third line, and if so, must be identical to the description on the first line.

Here is an example of one record from a FASTQ file:
```
@FSRRS4401BE7HA
tcagTTAAGATGGGAT
+
###EEEEEEEEE##E#
```

Where:
* `identifier` is `"FSRRS4401BE7HA"`
* `description` is also `"FSRRS4401BE7HA"`
* `sequence` is `"tcagTTAAGATGGGAT"`
* `quality` is `"###EEEEEEEEE##E#"`

## The `FASTQRecord`
`FASTQRecord`s optionally have the description repeated on the third line.
This can be toggled with `quality_header!(::Record, ::Bool)`:

```jldoctest qual
julia> record = parse(FASTQRecord, "@ILL01\nCCCGC\n+\nKM[^d");

julia> print(record)
@ILL01
CCCGC
+
KM[^d

julia> quality_header!(record, true); print(record)
@ILL01
CCCGC
+ILL01
KM[^d
```

```@docs
FASTQ.Record
```

## Qualities
Unlike `FASTARecord`s, a `FASTQRecord` contain quality scores, see the example above.

The quality string can be obtained using the `quality` method:
```jldoctest qual
julia> record = parse(FASTQRecord, "@ILL01\nCCCGC\n+\nKM[^d");

julia> quality(record)
"KM[^d"
```

Qualities are numerical values that are encoded by ASCII characters.
Unfortunately, multiple encoding schemes exist, although PHRED+33 is the most common.
The scores can be obtained using the `quality_scores` function, which returns an iterator of PHRED+33 scores:

```jldoctest qual
julia> collect(quality_scores(record))
5-element Vector{Int8}:
 42
 44
 58
 61
 67
```

If you want to decode the qualities using another scheme, you can use one of the predefined `QualityEncoding` objects.
For example, Illumina v 1.3 used PHRED+64:

```jldoctest qual
julia> collect(quality_scores(record, FASTQ.ILLUMINA13_QUAL_ENCODING))
5-element Vector{Int8}:
 11
 13
 27
 30
 36
```

Alternatively, `quality_scores` accept a name of the known quality encodings:

```jldoctest qual
julia> (collect(quality_scores(record, FASTQ.ILLUMINA13_QUAL_ENCODING)) ==
        collect(quality_scores(record, :illumina13)))
true
```

Lastly, you can create your own:

```@docs
QualityEncoding
```

### Reference:
```@docs
quality
quality_scores
quality_header!
```

## `FASTQReader` and `FASTQWriter`
`FASTQWriter` can optionally be passed the keyword `quality_header` to control whether or not to print the description on the third line (the one with `+`).
By default this is `nothing`, meaning that it will print the second header, if present in the record itself.

If set to a `Bool` value, the `Writer` will override the `Records`, without changing the records themselves.

### Reference:
```@docs
FASTQ
FASTQ.Reader
FASTQ.Writer
validate_fastq
```
