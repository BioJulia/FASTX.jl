```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# FASTA formatted files
__NB: First read the overview in the sidebar__

FASTA is a text-based file format for representing biological sequences.
A FASTA file stores a list of sequence records with name, description, and
sequence.

The template of a sequence record is:

```
>{description}
{sequence}
```

Where the "identifier" is the first part of the description up to the first whitespace
(or the entire description if there is no whitespace)

Here is an example of a chromosomal sequence:

```
>chrI chromosome 1
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC
CACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTA
```

## The `FASTARecord`
FASTA records are, by design, very lax in what they can contain.
They can contain almost arbitrary byte sequences, including invalid unicode, and trailing whitespace on their sequence lines, which will be interpreted as part of the sequence.
If you want to have more certainty about the format, you can either check the content of the sequences with a regex, or (preferably), convert them to the desired `BioSequence` type.

```@docs
FASTA.Record
```

### Reference:
```@docs
identifier
description
sequence
```

## `FASTAReader` and `FASTAWriter`
`FASTAWriter` can optionally be passed the keyword `width` to control the line width.
If this is zero or negative, it will write all record sequences on a single line.
Else, it will wrap lines to the given maximal width.

### Reference:
```@docs
FASTA.Reader
FASTA.Writer
```
