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

Here:
* The `identifier` is `"chrI"`
* The `description` is `"chrI chromosome 1"`, containing the identifier
* The sequence is the DNA sequence `"CCACA..."`

## The `FASTARecord`
FASTA records are, by design, very lax in what they can contain.
They can contain almost arbitrary byte sequences, including invalid Unicode, and trailing whitespace on their sequence lines, which will be interpreted as part of the sequence.
If you want to have more certainty about the format, you can either check the content of the sequences with a regex, or (preferably), convert them to the desired `BioSequence` type.

FASTX treats sequence positions as one-based byte positions and assumes ASCII
sequence symbols. In particular, `FASTA.Writer(width=...)` wraps after `width`
bytes—not Unicode characters—and may split a multi-byte UTF-8 symbol. For
multi-byte text, first obtain `sequence(String, record)` and perform
character-oriented work on a sequence constructed from that `String`.

```@docs
FASTA.Record
```

## `FASTAReader` and `FASTAWriter`
`FASTAWriter` can optionally be passed the keyword `width` to control the line width.
If this is zero or negative, it will write all record sequences on a single line.
Else, it will wrap lines to the given maximal width in bytes.

### Reference:
```@docs
FASTA
FASTA.Reader
FASTA.Writer
validate_fasta
```
