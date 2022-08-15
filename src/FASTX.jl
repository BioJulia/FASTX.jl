module FASTX

using StringViews: StringView
using BioSequences: BioSequence, LongSequence
using Automa: Automa

# using ReTest
#
# @testset "Inline test placeholder" begin
#     @test false
# end


"""
    identifier(record::Record)::AbstractString

Get the sequence identifier of `record`. The identifier is the description
before any whitespace. If the identifier is missing, return an empty string.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.

See also: [`description`](@ref), [`sequence`](@ref)

# Examples
```jldoctest
julia> record = parse(FASTA.Record, ">ident_here some descr \\nTAGA");

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
julia> record = parse(FASTA.Record, ">ident_here some descr \\nTAGA");

julia> description(record)
"ident_here some descr "
```
"""
function description end

"""
    sequence([::Type{S}], record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.BioSequence`, `AbstractString` or `String`.
If elided, `S` defaults to an `AbstractString` subtype.
If `part` argument is given, it returns the specified part of the sequence.

See also: [`identifier`](@ref), [`description`](@ref)

# Examples
```jldoctest
julia> record = parse(FASTQ.Record, "@read1\\nTAGA\\n+\\n;;]]");

julia> sequence(record)
"TAGA"

julia> sequence(LongDNA{2}, record)
4nt DNA Sequence:
TAGA
```
"""
function sequence end

# Get the length of the sequence part of a `FASTX` `Record` in number of bytes.
function seqlen end

const UTF8 = Union{AbstractVector{UInt8}, String, SubString{String}}

# line is nothing if the reader does not have line information after random IO access.
@noinline function throw_parser_error(
    data::Vector{UInt8},
    p::Integer,
    line::Union{Integer, Nothing}
)
    byte = data[p]
    # These bytes are printable in the Julia REPL as single chars e.g. "\t"
    bytestr = if byte in 0x07:0x13 || byte == 0x1b || byte in 0x20:0x7e
        ''' * Char(byte) * '''
    else
        repr(byte)
    end
    # These chars do not need escaping, e.g. '!', but not '\t'.
    bytestr = in(byte, 0x20:0x7e) ? bytestr : escape_string(bytestr)
    buf = IOBuffer()
    print(
        buf,
        "Error when parsing FASTX file. Saw unexpected byte ",
        bytestr
    )
    if line !== nothing
        print(buf, " on line ", string(line))

        # Compute column if possible, by looking at last '\n'.
        # it may not be possible because it may be past the data buffer `data`.
        lastnewline = findprev(isequal(UInt8('\n')), data, p)
        if lastnewline !== nothing
            col = p - lastnewline
            print(buf, " col ", string(col))
        end
    end
    error(String(take!(buf)))
end

# Truncate to at most `len` chars.
function truncate(s::AbstractString, len::Integer)
    if length(s) > len
        String(first(s, len-1) * '…')
    else
        return s
    end
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end

function appendfrom!(dst::Vector{UInt8}, dpos::Integer, src::Vector{UInt8}, spos::Integer, n::Integer)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copyto!(dst, dpos, src, spos, n)
    return dst
end

CONTEXT = Automa.CodeGenContext(
    generator=:goto,
    vars=Automa.Variables(:p, :p_end, :p_eof, :ts, :te, :cs, :data, :mem, :byte)
)

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

import .FASTA: FASTA, validate_fasta, Index, faidx, index!, extract, validate_fasta, seekrecord
import .FASTQ: FASTQ, quality, quality_scores, quality_header!, QualityEncoding, validate_fastq

function FASTA.Record(record::FASTQ.Record)
    slen = seqlen(record)
    dlen = record.description_len
    FASTA.Record(record.data[1:slen+dlen], record.identifier_len, dlen, slen)
end

"""
    copy!(::FASTA.Record, ::FASTQ.Record)

Copy the content of the FASTQ record into the FASTA record.
"""
function Base.copy!(dst::FASTA.Record, src::FASTQ.Record)
    dlen = src.description_len
    slen = seqlen(src)
    tlen = UInt(dlen + slen)
    dstdata = dst.data
    length(dstdata) < tlen && resize!(dstdata, tlen)
    copyto!(dstdata, 1, src.data, 1, tlen)
    dst.identifier_len = src.identifier_len
    dst.description_len = dlen
    dst.sequence_len = slen
    dst
end

Base.parse(::Type{T}, s::AbstractString) where {T <: Record} = parse(T, String(s))
Base.parse(::Type{T}, s::Union{String, SubString{String}}) where {T <: Record} = parse(T, codeunits(s))

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
    return S(sequence(record, part))
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
    validate_fasta,
    FASTQRecord,
    FASTQReader,
    FASTQWriter,
    validate_fastq,
    identifier,
    description,
    sequence,
    quality,
    quality_scores,
    quality_header!,
    QualityEncoding,
    Index,
    faidx,
    index!,
    extract,
    seekrecord

end # module
