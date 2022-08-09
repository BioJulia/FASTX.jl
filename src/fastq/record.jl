# FASTQ Record
# ============

"""
    FASTQ.Record

Mutable struct representing a FASTQ record as parsed from a FASTQ file.
The content of the record can be queried with the following functions:
`identifier`, `description`, `sequence`, `quality`
FASTQ records are un-typed, i.e. they are agnostic to what kind of data they contain.

See also: [`FASTQ.Reader`](@ref), [`FASTQ.Writer`](@ref)

# Examples
```jldoctest
julia> rec = parse(Record, "@ill r1\\nGGC\\n+\\njjk");

julia> identifier(rec)
"ill"

julia> description(rec)
"ill r1"

julia> sequence(rec)
"GGC"

julia> show(collect(quality(rec)))
[73, 73, 74]

julia> typeof(description(rec)) == typeof(sequence(rec)) isa AbstractString
true
```
"""
mutable struct Record
    # Contains: description, then sequence, then quality, then any noncoding bytes.
    # all rest, including newlines and the @ and + symbol, are not stored.
    # The second description after + must be identical to first description.
    # The quality is not corrected for offset, i.e. it is stored as it in the input file
    data::Vector{UInt8}
    
    identifier_len::Int32
    description_len::Int32

    # Top bit stores whether the description is repeated after the +
    has_description_seq_len::UInt
end

@inline seqlen(record::Record)::Int = (record.has_description_seq_len & (typemax(Int) % UInt)) % Int
has_extra_description(record::Record) = record.has_description_seq_len ≥ (typemin(Int) % UInt)

# Number of stored bytes in data field
filled(record::Record) = record.description_len + 2 * seqlen(record)

"""
    quality_header!(record::Record, x::Bool)

Set whether the record repeats its header on the quality comment line,
i.e. the line with `+`.

# Examples
```
julia> record = parse(FASTQ.Record, "@A B\\nT\\n+\\nJ");

julia> string(record)
"@A B\\nT\\n+\\nJ"

julia> quality_header!(record, true);

julia> string(record)
"@A B\\nT\\n+A B\\nJ"
```
"""
function quality_header!(record::Record, x::Bool)
    bits = record.has_description_seq_len
    y = if x
        bits | (typemin(Int) % UInt)
    else
        bits & (typemax(Int) % UInt)
    end
    record.has_description_seq_len = y
    record
end

"""
    FASTQ.Record()

Create the default FASTQ record.
"""
Record() = Record(UInt8[], 0, 0, 0)

function Base.empty!(record::Record)
    # Do not truncate the underlying data buffer
    record.identifier_len = 0
    record.description_len = 0
    record.has_description_seq_len = 0
    return record
end

function Base.parse(::Type{Record}, data::UTF8)
    record = Record()
    index!(record, data)
    return record
end

"""
    FASTQ.Record(description, sequence, quality; offset=33)

Create a FASTQ record from `description`, `sequence` and `quality`.
Arguments:
* `description::AbstractString`
* `sequence::Union{AbstractString, BioSequence}`,
* `quality::Union{AbstractString, Vector{<:Number}}`
* Keyword argument `offset` (if `quality isa Vector`): PHRED offset
"""
function Record(
    description::AbstractString,
    sequence::Union{AbstractString, BioSequence},
    quality::AbstractString
)
    seqlen = sequence isa AbstractString ? ncodeunits(sequence) : length(sequence)
    if seqlen != ncodeunits(quality)
        throw(ArgumentError("Byte length of sequence doesn't match codeunits of quality"))
    end
    buf = IOBuffer()
    print(buf,
        '@', description, '\n',
        sequence, "\n+\n",
        quality
    )
    parse(Record, take!(buf))
end

function Record(
    description::AbstractString,
    sequence::Union{AbstractString, BioSequence},
    quality::Vector{<:Number};
    offset::Integer=33
)
    ascii_quality = String([UInt8(q + offset) for q in quality])
    Record(description, sequence, ascii_quality)
end

function Base.:(==)(record1::Record, record2::Record)
    record1.description_len == record2.description_len || return false
    filled1 = filled(record1)
    filled1 == filled(record2) || return false
    (data1, data2) = (record1.data, record2.data)
    GC.@preserve data1 data2 begin
        return memcmp(pointer(data1), pointer(data2), filled1) == 0
    end
end

function Base.copy(record::Record)
    return Record(
        record.data[1:filled(record)],
        record.identifier_len,
        record.description_len,
       record.has_description_seq_len
    )
end

function Base.write(io::IO, record::Record)
    data = record.data
    len = UInt(seqlen(record))
    desclen = UInt(record.description_len)
    GC.@preserve data begin
        # Header line
        nbytes = write(io, UInt8('@'))
        nbytes += unsafe_write(io, pointer(data), desclen)
        nbytes += write(io, '\n')
        # Sequence
        nbytes += unsafe_write(io, pointer(data) + desclen, len)
        # + line, with additional description if applicable
        nbytes += write(io, UInt8('\n'), UInt8('+'))
        if has_extra_description(record)
            nbytes += unsafe_write(io, pointer(data), desclen)
        end
        # Quality
        nbytes += write(io, '\n')
        nbytes += unsafe_write(io, pointer(data) + desclen + len, len)
    end
    return nbytes
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    println(io, "FASTQ.Record:")
    println(io, "  description: ", description(record))
    println(io, "     sequence: ", sequence(String, record))
    print(io,   "      quality: ", quality(String, record))
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end


# Accessor functions
# ------------------

function quality_indices(record::Record, part::UnitRange{<:Integer})
    start, stop = first(part), last(part)
    (start < 1 || stop > seqlen(record)) && throw(BoundsError(record, start:stop))
    offset = record.description_len + seqlen(record)
    start+offset:stop+offset
end

"""
    quality([T::Type{String, StringView}], record::FASTQ.Record, [part::UnitRange])

Get the ASCII quality of `record` at positions `part` as type `T`.
If not passed, `T` defaults to `StringView`.
If not passed, `part` defaults to the entire quality string.

# Examples
```jldoctest
julia> rec = parse(FASTQ.Record, "@hdr\\nUAGUCU\\n+\\nCCDFFG");

julia> quality(rec)
"CCDFFG"

julia> quality isa AbstractString
true
```
"""
function quality(record::Record, part::UnitRange{<:Integer}=1:seqlen(record))
    quality(StringView, record, part)
end

function quality(::Type{String}, record::Record, part::UnitRange{<:Integer}=1:seqlen(record))
    String(record.data[quality_indices(record, part)])
end

function quality(::Type{StringView}, record::Record, part::UnitRange{<:Integer}=1:seqlen(record))
    StringView(view(record.data, quality_indices(record, part)))
end

function quality_scores(record::Record, part::UnitRange{<:Integer}=1:seqlen(record))
    quality_scores(record, DEFAULT_ENCODING, part)
end

"""
    quality_scores(record::FASTQ.Record, [encoding::QualityEncoding], [part::UnitRange])

Get an iterator of PHRED base quality scores of `record` at positions `part`.
This iterator is corrupted if the record is mutated.
By default, `part` is the whole sequence.
By default, the encoding is PHRED33 Sanger encoding, but may be specified with a `QualityEncoding` object
"""
function quality_scores(record::Record, encoding::QualityEncoding, part::UnitRange{<:Integer}=1:seqlen(record))
    start, stop = first(part), last(part)
    (start < 1 || stop > seqlen(record)) && throw(BoundsError(record, start:stop))
    data = record.data
    offset = record.description_len + seqlen(record) 
    return Iterators.map(offset+start:offset+stop) do i
        v = data[i]
        decode_quality(encoding, v)
    end
end

"""
    quality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}

Get an iterator of base quality of the slice `part` of `record`'s quality.

The `encoding_name` can be either `:sanger`, `:solexa`, `:illumina13`, `:illumina15`, or `:illumina18`.
"""
function quality_scores(
    record::Record,
    encoding_name::Symbol,
    part::UnitRange{<:Integer}=1:seqlen(record)
)
    encoding = (
        encoding_name == :sanger     ?     SANGER_QUAL_ENCODING :
        encoding_name == :solexa     ?     SOLEXA_QUAL_ENCODING :
        encoding_name == :illumina13 ? ILLUMINA13_QUAL_ENCODING :
        encoding_name == :illumina15 ? ILLUMINA15_QUAL_ENCODING :
        encoding_name == :illumina18 ? ILLUMINA18_QUAL_ENCODING :
        throw(ArgumentError("quality encoding ':$(encoding_name)' is not supported")))
    quality_scores(record, encoding, part)
end

function Base.hash(record::Record, h::UInt)
    h = hash(record.description_len, h)
    hash(view(record.data, filled(record)), h)
end
