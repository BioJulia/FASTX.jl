# FASTQ Read
# ==========
#
# Design discussions in https://github.com/BioJulia/FASTX.jl/pull/35

struct FASTQRead{A <: BioSequences.Alphabet}
    identifier::String
    description::String
    sequence::BioSequences.LongSequence{A}
    quality::Vector{UInt8} # in raw PHRED scores, not offset
end

"""
    FASTQ.Read(record::FASTQ.Record, offset::Integer=33)
    FASTQ.Read(record::FASTQ.Record, encoding_name::Symbol)

Create a FASTQ read object from a FASTQ.Record.
The FASTQ.Read has fields: 
    identifier::String
    description::String
    sequence::BioSequences.LongSequence
    quality::Vector{UInt8} # in raw PHRED scores, not offset
"""
function Read(record::FASTQ.Record, offset::Integer=33)
    FASTQRead(
        identifier(record),
        description(record),
        sequence(record),
        quality(record, offset)
    )
end

function Read(record::FASTQ.Record, encoding_name::Symbol)
    FASTQRead(
        identifier(record),
        description(record),
        sequence(record),
        quality(record, encoding_name)
    )
end

"""
    sequence(read::FASTQ.FASTQRead)
Get the sequence of a FASTQ read. Same as read.sequence
"""
function sequence(read::FASTQ.FASTQRead)
    read.sequence
end

"""
    quality(read::FASTQ.FASTQRead)
Get the quality of a FASTQ read (Vector{UInt8}). Same as read.quality
"""
function quality(read::FASTQ.FASTQRead)
    read.quality
end

"""
    length(read::FASTQ.FASTQRead)
Get the length of a FASTQ read.
"""
function Base.length(read::FASTQ.FASTQRead)
    length(read.sequence)
end

"""
     Base.getindex(record::Record, i::UnitRange{Int})
Subset a FASTQRead using string syntax. Eg. read[3:7] or read[3]
"""
function Base.getindex(read::FASTQ.FASTQRead, i::UnitRange{<:Integer})
    return FASTQRead(
        read.identifier,
        read.description,
        read.sequence[i],
        read.quality[i]
    )
end

function Base.getindex(read::FASTQRead, i::Integer)
    read[i:i]
end

function Base.show(io::IO, read::FASTQ.FASTQRead)
    print(io, summary(read), ':')
    println(io)
    println(io, "   identifier: ", read.identifier)
    println(io, "  description: ", read.description)
    println(io, "     sequence: ", String(read.sequence))
    print(io, "      quality: ", Int64.(quality(read)))
end
