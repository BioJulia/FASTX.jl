```@meta
CurrentModule = FASTX
DocTestSetup = quote
    using FASTX
end
```

# Public API Reference

## Contents

```@contents
Pages = ["public.md"]
```

```@index
Pages = ["public.md"]
```

## FASTA API

The following methods and types are provided by the FASTA
submodule for public use. They are not exported as in general
using FASTX requires qualifying the submodule (FASTA or FASTQ)
that you are using.

```@docs
FASTA.Reader
FASTA.Writer
FASTA.Record
FASTA.hasidentifier
FASTA.identifier
FASTA.hasdescription
FASTA.description
FASTA.hassequence
FASTA.sequence
FASTA.seqlen
```

## FASTQ API

The following methods and types are provided by the FASTQ
submodule for public use. They are not exported as in general
using FASTX requires qualifying the submodule (FASTA or FASTQ)
that you are using.

```@docs
FASTQ.Reader
FASTQ.Writer
FASTQ.Record
FASTQ.hasidentifier
FASTQ.identifier
FASTQ.hasdescription
FASTQ.description
FASTQ.hassequence
FASTQ.sequence
FASTQ.seqlen
FASTQ.hasquality
FASTQ.quality
FASTQ.SANGER_QUAL_ENCODING
FASTQ.SOLEXA_QUAL_ENCODING
FASTQ.ILLUMINA13_QUAL_ENCODING
FASTQ.ILLUMINA15_QUAL_ENCODING
FASTQ.ILLUMINA18_QUAL_ENCODING
```