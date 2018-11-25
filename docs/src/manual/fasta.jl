# # IO - FASTA formatted files
#
# FASTA is a text-based file format for representing biological sequences.
# A FASTA file stores a list of sequence records with name, description, and
# sequence.
#
# The template of a sequence record is:
#
# ```
# >{name} {description}?
# {sequence}
# ```
#
# Here is an example of a chromosomal sequence:
#
# ```
# >chrI chromosome 1
# CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC
# CACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG
# ```
#
# ## Readers and Writers
# The reader and writer for FASTA formatted files, are found within the
# `BioSequences.FASTA` submodule.
#
#md # ```@docs
#md # FASTA.Reader
#md # FASTA.Writer
#md # ```
#
# They can be created with IOStreams:

using FASTX

r = FASTA.Reader(open("MyInput.fasta", "r"))
w = FASTA.Writer(open("MyFile.fasta", "w"))

# Alternatively, `Base.open` is overloaded with a method for conveinience:

r = open(FASTA.Reader, "MyInput.fasta")
w = open(FASTA.Writer, "MyFile.fasta")

# Usually sequence records will be read sequentially from a file by iteration.

reader = open(FASTA.Reader, "hg38.fa")
for record in reader
    ## Do something
end
close(reader)

# But if the FASTA file has an auxiliary index file formatted in fai, the reader
# supports random access to FASTA records, which would be useful when accessing
# specific parts of a huge genome sequence:

reader = open(FASTAReader, "sacCer.fa", index="sacCer.fa.fai")
chrIV = reader["chrIV"]  # directly read sequences called chrIV.
close(reader)

# Reading in a sequence from a FASTA formatted file will give you a variable of
# type `FASTA.Record`.
#
#md # ```@docs
#md # FASTA.Record
#md # ```
#
# Various getters and setters are available for `FASTA.Record`s:
#
#md # ```@docs
#md # FASTA.hasidentifier
#md # FASTA.identifier
#md # FASTA.hasdescription
#md # FASTA.description
#md # FASTA.hassequence
#md # FASTA.sequence(record::FASTA.Record, [part::UnitRange{Int}])
#md # ```
#
# To write a `BioSequence` to FASTA file, you first have to create a `FASTA.Record`:

using BioSequences
x = dna"aaaaatttttcccccggggg"
rec = FASTA.Record("MySeq", x)
w = open(FASTA.Writer, "MyFile.fasta")
write(w, rec)
close(w)

# As always with julia IO types, remember to close your file readers and writer
# after you are finished.

# Using `open` with a do-block can help ensure you close a stream after you are
# finished.

open(FASTA.Reader, "hg38.fa") do reader
    for record in reader
        ## Do something
    end
end

