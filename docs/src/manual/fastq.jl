# # IO - FASTQ formatted files
#
# FASTQ is a text-based file format for representing DNA sequences along with
# qualities for each base.
# A FASTQ file stores a list of sequence records in the following format:
#
# ```
# @{name} {description}?
# {sequence}
# +
# {qualities}
# ```
#
# Here is an example of one record from a FASTQ file:
#
# ```
# @FSRRS4401BE7HA
# tcagTTAAGATGGGAT
# +
# ###EEEEEEEEE##E#
# ```
#
# ## Readers and Writers
# The reader and writer for FASTQ formatted files, are found within the
# `BioSequences.FASTQ` module.
#
#md ```@docs
#md FASTQ.Reader
#md FASTQ.Writer
#md ```
#
# They can be created with IOStreams:

r = FASTQ.Reader(open("MyInput.fastq", "r"))
w = FASTQ.Writer(open("MyFile.fastq", "w"))

# Note that `FASTQ.Reader` does not support line-wraps within sequence and quality.
# Usually sequence records will be read sequentially from a file by iteration.

using FASTX
reader = FASTQ.Reader(open("hg38.fastq", "r"))
for record in reader
    ## Do something
end
close(reader)

# Reading in a record from a FASTQ formatted file will give you a variable of
# type `FASTQ.Record`.

#md ```@docs
#md FASTQ.Record
#md ```

# Various getters and setters are available for `FASTQ.Record`s:
#
#md ```@docs
#md FASTQ.hasidentifier
#md FASTQ.identifier
#md FASTQ.hasdescription
#md FASTQ.description
#md FASTQ.hassequence
#md FASTQ.sequence(record::FASTQ.Record, [part::UnitRange{Int}])
#md FASTQ.hasquality
#md FASTQ.quality
#md ```
#
# To write a `BioSequence` to FASTQ file, you first have to create a `FASTQ.Record`:
#
#md ```@docs
#md FASTQ.Record(identifier::AbstractString, description::Union{AbstractString,Nothing}, sequence, quality::Vector; offset=33)
#md ```
#
# As always with julia IO types, remember to close your file readers and writer
# after you are finished.
#
# Using `open` with a do-block can help ensure you close a stream after you are
# finished.

open(FASTQ.Reader, "hg38.fastq") do reader
    for record in reader
        ## Do something
    end
end