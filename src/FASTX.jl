module FASTX

export
    FASTA,
    FASTQ,
    identifier,
    description,
    header,
    sequence,
    quality,
    quality_iter,
    transcribe

using StringViews: StringView
using BioSequences: BioSequence, LongSequence

# Generic methods

"""
    identifier(record::Record)::StringView

Get the sequence identifier of `record`. The identifier is the header
before any whitespace. If the identifier is missing, return an empty string.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.
"""
function identifier end

"""
    description(record::Record)::StringView

Get the description of `record`. The description is the entire header line.
Returns an `AbstractString` view into the record.
"""
function description end

"""
    sequence([::Type{S}], record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.BioSequence` or `String`.
If elided, `S` defaults to `StringView`.
If `part` argument is given, it returns the specified part of the sequence.

!!! note
    This method makes a new sequence object every time.
    If you have a sequence already and want to fill it with the sequence
    data contained in a record, you can use `Base.copyto!`.
"""
function sequence end
function seqlen end

const UTF8 = Union{AbstractVector{UInt8}, String, SubString{String}}

const header = description

include("fasta/fasta.jl")
include("fastq/fastq.jl")

const Record = Union{FASTA.Record, FASTQ.Record}

"Get the indices of `data` that correspond to sequence indices `part`"
function seq_data_part(record::Record, part::AbstractUnitRange)
    start, stop = first(part), last(part)
    (start < 1 || stop > seqlen(record)) && throw(BoundsError(record, start:stop))
    Int(record.description_len) + start:Int(record.description_len) + stop
end

sequence(record::Record, part::UnitRange{Int}=1:seqlen(record)) = sequence(StringView, record, part)

function sequence(::Type{StringView}, record::Record, part::UnitRange{Int}=1:seqlen(record))
    return StringView(view(record.data, seq_data_part(record, part)))
end

function sequence(
    ::Type{S},
    record::Record,
    part::UnitRange{Int}=1:seqlen(record)
)::S where S <: BioSequence
    return S(@view(record.data[seq_data_part(record, part)]))
end

function sequence(
    ::Type{String},
    record::Record,
    part::UnitRange{Int}=1:seqlen(record)
)::String
    return String(record.data[seq_data_part(record, part)])
end

function Base.copy!(dest::LongSequence, src::Record)
    resize!(dest, UInt(seqlen(src)))
    copyto!(dest, 1, src, 1, seqlen(src))
end

"""
    Base.copyto!(dest::BioSequences.LongSequence, src::Record)

Copy all of the sequence data from the fastq record `src` to a biological
sequence `dest`. `dest` must have a length greater or equal to the length of
the sequence represented in the fastq record. The first n elements of `dest` are
overwritten, the other elements are left untouched.
"""
function Base.copyto!(dest::LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, seqlen(src))
end

import .FASTA
import .FASTQ
import .FASTQ: quality, hasquality, quality_iter

# TODO: Check this
#=
function FASTA.Record(record::FASTQ.Record)
    FASTQ.checkfilled(record)
    slice = (first(record.identifier) - 1):last(record.sequence)
    newdata = @inbounds record.data[slice]
    newdata[1] = 0x3E
    FASTA.Record(newdata)
end

"""
    transcribe(in::FASTQ.Reader, out::FASTA.Writer)

Convert a FASTQ file to a FASTA file.
"""
function transcribe(in::FASTQ.Reader, out::FASTA.Writer)
    buff_record = FASTQ.Record()
    while !eof(in)
        read!(in, buff_record)
        write(out, FASTA.Record(buff_record))
    end
end
=#

end # module
