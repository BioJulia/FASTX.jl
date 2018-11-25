# FASTQ

## FASTQ Reader

```@docs
FASTQ.Reader
```

## FASTQ Record

```@docs
FASTQ.Record
FASTQ.identifier
FASTQ.hasidentifier
FASTQ.description
FASTQ.hasdescription
FASTQ.sequence
FASTQ.hassequence
FASTQ.quality
FASTQ.hasquality
FASTQ.QualityEncoding
FASTQ.SANGER_QUAL_ENCODING
FASTQ.SOLEXA_QUAL_ENCODING
FASTQ.ILLUMINA13_QUAL_ENCODING
FASTQ.ILLUMINA15_QUAL_ENCODING
FASTQ.ILLUMINA18_QUAL_ENCODING 
```

### BioCore methods

The following BioCore methods will work with
`FASTQ.Record` types.

```@docs
BioCore.isfilled
BioCore.seqname
BioCore.hasseqname
BioCore.sequence
BioCore.hassequence
``` 

## FASTQ Writer

```@docs
FASTQ.Writer
```

