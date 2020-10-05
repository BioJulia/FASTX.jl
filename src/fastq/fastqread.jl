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

Create a FASTQ read object from a FASTQ.Record.
The FASTQ.Read has fields: 
    identifier::String
    description::String
    sequence::BioSequences.LongDNASeq
    quality::Vector{UInt8} # in raw PHRED scores, not offset
"""
function Read(record::FASTQ.Record, offset::Integer=33)
    FASTQRead(
        identifier(record),
        description(record),
        sequence(BioSequences.LongDNASeq, record),
        quality(record, offset)
    )
end

function sequence(read::FASTQ.FASTQRead)
    read.sequence
end

function quality(read::FASTQ.FASTQRead)
    read.quality
end

function Base.length(read::FASTQ.FASTQRead)
    length(read.sequence)
end

function Base.getindex(read::FASTQ.FASTQRead, i::UnitRange{Int})
    return FASTQRead(
        read.identifier,
        read.description,
        read.sequence[i],
        read.quality[i]
    )
end

function Base.getindex(read::FASTQRead, i::Int)
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
