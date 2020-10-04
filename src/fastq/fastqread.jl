# FASTQ Read
# ==========
#
# Design discussions in https://github.com/BioJulia/FASTX.jl/pull/35

mutable struct FASTQRead
    identifier::String
    description::String
    sequence::BioSequences.LongDNASeq
    quality::Vector{UInt8} # in raw PHRED scores, not offset
end


function Read(record::Record, offset::Integer=33)
    FASTQRead(
        identifier(record),
        description(record),
        sequence(BioSequences.LongDNASeq, record),
        quality(record, offset)
    )
end

function sequence(read::Read)
    read.sequence
end

function quality(read::Read)
    read.quality
end

function Base.length(read)
    length(read.lequence)
end

function Base.getindex(read::Read, i::UnitRange{Int})
    return Read(
        read.identifier,
        read.description,
        read.sequence[i],
        read.quality[i]
    )
end

function Base.getindex(read::Read, i::Int)
    read[i:i]
end
