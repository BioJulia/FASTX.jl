module FASTX

using StringViews: StringView
using BioSequences: BioSequence, LongSequence

"""
    identifier(record::Record)::AbstractString

Get the sequence identifier of `record`. The identifier is the description
before any whitespace. If the identifier is missing, return an empty string.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.

See also: [`description`](@ref), [`sequence`](@ref)

# Examples
```jldoctest
julia> record = parse(FASTA.Record, ">ident_here some descr \nTAGA");

julia> identifier(record)
"ident_here"
```
"""
function identifier end

"""
    description(record::Record)::AbstractString

Get the description of `record`. The description is the entire header line, minus the
leading `>` or `@` symbols for FASTA/FASTQ records, respectively, including trailing whitespace.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.

See also: [`identifier`](@ref), [`sequence`](@ref)

# Examples
```jldoctest
julia> record = parse(FASTA.Record, ">ident_here some descr \nTAGA");

julia> description(record)
"some descr "
```
"""
function description end

"""
    sequence([::Type{S}], record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.BioSequence`, `AbstractString` or `String`.
If elided, `S` defaults to an `AbstractString`.
If `part` argument is given, it returns the specified part of the sequence.

See also: [`identifier`](@ref), [`description`](@ref)

# Examples
```jldoctest
julia> record = parse(FASTQ.Record, "@read1\nTAGA\n+\n;;]]");

julia> sequence(record)
"TAGA"

julia> sequence(LongDNA{2}, record)
4nt DNA Sequence:
TAGA
```
"""
function sequence end

"""
    seqlen(::Record)::Int

Get the length of the sequence part of a `FASTX` `Record` in number of bytes.
"""
function seqlen end

const UTF8 = Union{AbstractVector{UInt8}, String, SubString{String}}

include("fasta/fasta.jl")
include("fastq/fastq.jl")

const Record = Union{FASTA.Record, FASTQ.Record}

# Generic methods

function identifier(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.identifier_len)))
end

function description(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.description_len)))
end

import .FASTA: FASTA, validate_fasta, Index, faidx, extract
import .FASTQ: FASTQ, quality, quality_scores, quality_header!, QualityEncoding

function FASTA.Record(record::FASTQ.Record)
    ilen = record.identifier_len
    dlen = record.description_len
    slen = seqlen(record)
    tlen = UInt(dlen + slen)
    FASTA.Record(record.data[1:tlen], ilen, dlen, slen)
end

"""
    FASTA.Record!(::FASTQ.Record)

Convert the `FASTQ.Record` to a `FASTA.Record`, taking control of the underlying
data. The FASTQ record cannot be used after this operation.
"""
function FASTA.Record!(record::FASTQ.Record)
    FASTA.Record(record.data, record.identifier_len, record.description_len, seqlen(record))
end

Base.parse(::Type{T}, s::AbstractString) where {T <: Record} = parse(T, String(s))

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
    return S(sequence(record))
end

# Special method for LongSequence: Can operate on bytes directly
# and more efficiently
function sequence(
    ::Type{S},
    record::Record,
    part::UnitRange{Int}=1:seqlen(record)
)::S where S <: LongSequence
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

function Base.copyto!(dest::LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, seqlen(src))
end

function Base.copyto!(dest::LongSequence, doff, src::Record, soff, N)
    # This check is here to prevent boundserror when indexing src.sequence
    iszero(N) && return dest
    return copyto!(dest, doff, src.data, Int(src.description_len) + soff, N)
end

const FASTARecord = FASTA.Record
const FASTQRecord = FASTQ.Record
const FASTAReader = FASTA.Reader
const FASTQReader = FASTQ.Reader
const FASTAWriter = FASTA.Writer
const FASTQWriter = FASTQ.Writer

export
    FASTA,
    FASTQ,
    FASTARecord,
    FASTAReader,
    FASTAWriter,
    FASTQRecord,
    FASTQReader,
    FASTQWriter,
    identifier,
    description,
    sequence,
    quality,
    quality_scores,
    quality_header!,
    QualityEncoding,
    seqlen,
    Index,
    faidx,
    extract

end # module
