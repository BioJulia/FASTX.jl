# FASTQ Record
# ============

mutable struct Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    identifier::UnitRange{Int}
    description::UnitRange{Int}
    sequence::UnitRange{Int}
    quality::UnitRange{Int}
end

"""
    FASTQ.Record()

Create an unfilled FASTQ record.
"""
function Record()
    return Record(UInt8[], 1:0, 1:0, 1:0, 1:0, 1:0)
end

"""
    FASTQ.Record(data::Vector{UInt8})

Create a FASTQ record object from `data`.

This function verifies and indexes fields for accessors.

!!! warning
    Note that the ownership of `data` is transferred to a new record object.
    Editing the input data will edit the record, and is not advised after 
    construction of the record.
"""
function Record(data::Vector{UInt8})
    record = Record(data, 1:0, 1:0, 1:0, 1:0, 1:0)
    index!(record)
    return record
end

"""
    FASTQ.Record(str::AbstractString)

Create a FASTQ record object from `str`.

This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return Record(Vector{UInt8}(str))
end

Base.parse(::Record, str::AbstractString) = Record(str)

"""
    FASTQ.Record(identifier, sequence, quality; offset=33)

Create a FASTQ record from `identifier`, `sequence` and `quality`.
"""
function Record(identifier::AbstractString, sequence, quality::Vector; offset=33)
    return Record(identifier, nothing, sequence, quality, offset=offset)
end

"""
    FASTQ.Record(identifier, description, sequence, quality; offset=33)

Create a FASTQ record from `identifier`, `description`, `sequence` and `quality`.
"""
function Record(identifier::AbstractString, description::Union{AbstractString,Nothing}, sequence, quality::Vector; offset=33)
    if length(sequence) != length(quality)
        throw(ArgumentError("the length of sequence doesn't match the length of quality"))
    end
    buf = IOBuffer()
    print(buf, '@', identifier)
    if description != nothing
        print(buf, ' ', description)
    end
    print(buf, '\n')
    print(buf, sequence, '\n')
    print(buf, "+\n")
    ascii_quality = convert(Vector{UInt8}, quality .+ offset)
    write(buf, ascii_quality, '\n')
    return Record(take!(buf))
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.identifier,
        record.description,
        record.sequence,
        record.quality)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   identifier: ", hasidentifier(record) ? identifier(record) : "<missing>")
        println(io, "  description: ", hasdescription(record) ? description(record) : "<missing>")
        println(io, "     sequence: ", hassequence(record) ? sequence(String, record) : "<missing>")
          print(io, "      quality: ", hasquality(record) ? quality(record) : "<missing>")
    else
        print(io, " <not filled>")
    end
end

function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    record.quality = 1:0
    return record
end

function BioGenerics.isfilled(record::Record)
    return !isempty(record.filled)
end


# Accessor functions
# ------------------

"""
    identifier(record::Record)::Union{String,Nothing}

Get the sequence identifier of `record`.

!!! note
    Returns `nothing` if the record has no identifier.
"""
function identifier(record::Record)::Union{String,Nothing}
    checkfilled(record)
    if !hasidentifier(record)
        return nothing
    end
    return String(record.data[record.identifier])
end

"""
    hasidentifier(record::Record)

Checks whether or not the `record` has an identifier.
"""
function hasidentifier(record::Record)
    return isfilled(record)
end


"""
    description(record::Record)::Union{String, Nothing}

Get the description of `record`.

!!! note
    Returns `nothing` if `record` has no description.
"""
function description(record::Record)::Union{String, Nothing}
    checkfilled(record)
    if !hasdescription(record)
        nothing
    end
    return String(record.data[record.description])
end

"""
    hasdescription(record::Record)

Checks whether or not the `record` has a description.
"""
function hasdescription(record)
    return isfilled(record) && record.description != 1:0
end

"""
    sequence(::Type{S}, record::Record, [part::UnitRange{Int}])

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.BioSequence` or `String`.
If `part` argument is given, it returns the specified part of the sequence.
"""
function sequence(::Type{S}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::S where S <: BioSequences.BioSequence
    checkfilled(record)
    seqpart = record.sequence[part]
    return S(record.data, first(seqpart), last(seqpart))
end

"""
    sequence(::Type{String}, record::Record, [part::UnitRange{Int}])::String

Get the sequence of `record` as a String.
If `part` argument is given, it returns the specified part of the sequence.
"""
function sequence(::Type{String}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::String
    checkfilled(record)
    return String(record.data[record.sequence[part]])
end

"""
    sequence(record::Record, [part::UnitRange{Int}])::BioSequences.DNASequence
Get the sequence of `record`.
"""
function sequence(record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::BioSequences.DNASequence
    return sequence(BioSequences.DNASequence, record, part)
end

"""
    hassequence(record::Record)

Checks whether or not a sequence record contains a sequence.

!!! note
    Zero-length sequences are allowed in records.
"""
function hassequence(record::Record)
    # zero-length sequence may exist
    return isfilled(record)
end

"""
    quality(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}

Get the base quality of `record`.
"""
function quality(record::Record, offset::Integer=33, part::UnitRange{Int}=1:lastindex(record.quality))::Vector{UInt8}
    checkfilled(record)
    quality = record.data[record.quality[part]]
    for i in 1:lastindex(part)
        # TODO: Checked arithmetic?
        @inbounds quality[i] -= offset
    end
    return quality
end

"""
    quality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}

Get the base quality of `record` by decoding with `encoding_name`.

The `encoding_name` can be either `:sanger`, `:solexa`, `:illumina13`, `:illumina15`, or `:illumina18`.

!!! note
    Returns `nothing` if the record has no quality string.
"""
function quality(record::Record, encoding_name::Symbol, part::UnitRange{Int}=1:lastindex(record.quality))::Vector{UInt8}
    checkfilled(record)
    encoding = (
        encoding_name == :sanger     ?     SANGER_QUAL_ENCODING :
        encoding_name == :solexa     ?     SOLEXA_QUAL_ENCODING :
        encoding_name == :illumina13 ? ILLUMINA13_QUAL_ENCODING :
        encoding_name == :illumina15 ? ILLUMINA15_QUAL_ENCODING :
        encoding_name == :illumina18 ? ILLUMINA18_QUAL_ENCODING :
        throw(ArgumentError("quality encoding ':$(encoding_name)' is not supported")))
    quality = Vector{UInt8}(undef, length(part))
    if !isempty(part)
        qpart = record.quality[part]
        check_quality_string(encoding, record.data, first(qpart), last(qpart))
        decode_quality_string!(encoding, record.data, quality, first(qpart), last(qpart))
    end
    return quality
end

"""
    hasquality(record::Record)

Check whether the given FASTQ `record` has a quality string.
"""
function hasquality(record::Record)
    return isfilled(record)
end

function BioGenerics.seqname(record::Record)
    return identifier(record)
end

function BioGenerics.hasseqname(record::Record)
    return hasidentifier(record)
end

function BioGenerics.sequence(record::Record)
    return sequence(record)
end

function BioGenerics.sequence(::Type{S}, record::Record) where S <: BioSequences.BioSequence
    return sequence(S, record)
end

function BioGenerics.hassequence(record::Record)
    return hassequence(record)
end


# Helper functions
# ----------------

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled FASTQ record"))
    end
end
