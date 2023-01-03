module BioSequencesExt

import FASTX: FASTARecord, FASTQRecord, sequence, Record, seqsize, seq_data_part
using BioSequences: BioSequence, LongSequence

function sequence(
    ::Type{S},
    record::Record,
    part::UnitRange{Int}=1:seqsize(record)
)::S where S <: BioSequence
    return S(sequence(record, part))
end

# Special method for LongSequence: Can operate on bytes directly
# and more efficiently
function sequence(
    ::Type{S},
    record::Record,
    part::UnitRange{Int}=1:seqsize(record)
)::S where S <: LongSequence
    return S(@view(record.data[seq_data_part(record, part)]))
end

function Base.copy!(dest::LongSequence, src::Record)
    resize!(dest, UInt(seqsize(src)))
    copyto!(dest, 1, src, 1, seqsize(src))
end

function Base.copyto!(dest::LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, seqsize(src))
end

function Base.copyto!(dest::LongSequence, doff, src::Record, soff, N)
    # This check is here to prevent boundserror when indexing src.sequence
    iszero(N) && return dest
    return copyto!(dest, doff, src.data, Int(src.description_len) + soff, N)
end

function FASTARecord(description::AbstractString, sequence::BioSequence)
    buf = IOBuffer()
    print(buf, '>', description, '\n')
    # If the sequence is empty, we need to print a newline in order to not
    # have the FASTA file truncated, thus invalid
    print(buf, isempty(sequence) ? '\n' : sequence)
    return parse(FASTARecord, take!(buf))
end

function FASTQRecord(
    description::AbstractString,
    sequence::BioSequence,
    quality::AbstractString
)
    if length(sequence) != ncodeunits(quality)
        throw(ArgumentError("Byte length of sequence doesn't match codeunits of quality"))
    end
    buf = IOBuffer()
    print(buf,
        '@', description, '\n',
        sequence, "\n+\n",
        quality
    )
    parse(FASTQRecord, take!(buf))
end

function FASTQRecord(
    description::AbstractString,
    sequence::BioSequence,
    quality::Vector{<:Number};
    offset::Integer=33
)
    ascii_quality = String([UInt8(q + offset) for q in quality])
    FASTQRecord(description, sequence, ascii_quality)
end

end
