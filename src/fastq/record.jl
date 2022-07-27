# FASTQ Record
# ============

mutable struct Record
    # Contains: description, then sequence, then quality, then any noncoding bytes.
    # all rest, including newlines and the @ and + symbol, are not stored.
    # The second description after + must be identical to first description
    data::Vector{UInt8}
    
    identifier_len::Int32
    description_len::Int32

    # Top bit stores whether the description is repeated after the +
    has_description_seq_len::UInt
end

@inline seqlen(record::Record) = (record.has_description_seq_len & (typemax(Int) % UInt)) % Int
has_extra_description(record::Record) = record.has_description_seq_len ≥ (typemin(Int) % UInt)

# Number of stored bytes in data field
filled(record::Record) = record.description_len + 2 * seqlen(record)

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


"""
    FASTQ.Record(data::Union{Vector{UInt8, String, SubString{String}}})

Create a FASTQ record object from `data`.

This function verifies and indexes fields for accessors.
Note that this function allocates a new array.
To parse a record in-place, use `index!(record, data)`
"""
Record(data::AbstractString) = Record(String(data))
Base.parse(::Type{Record}, s::AbstractString) = Record(s)

function Record(data::UTF8)
    record = Record()
    index!(record, data)
    return record
end

# TODO: This could be done more efficiently.
"""
    FASTQ.Record(description, sequence, quality; offset=33)

Create a FASTQ record from `description`, `sequence` and `quality`.
"""
function Record(description::AbstractString, sequence, quality::Vector{<:Number}; offset::Integer=33)
    if length(sequence) != length(quality)
        throw(ArgumentError("the length of sequence doesn't match the length of quality"))
    end
    buf = IOBuffer()
    print(buf, '@', description, '\n')
    print(buf, sequence, "\n+\n")
    ascii_quality = [UInt8(q + offset) for q in quality]
    write(buf, ascii_quality, '\n')
    return Record(take!(buf))
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
    seqlen = UInt(seqlen(record))
    desclen = UInt(sequence.description_len)
    GC.@preserve data begin
        # Header line
        nbytes = write(io, UInt8('@'))
        nbytes += unsafe_write(io, pointer(data), desclen)
        nbytes += write(io, '\n')
        # Sequence
        nbytes += unsafe_write(io, pointer(data) + desclen, seqlen)
        # + line, with additional description if applicable
        nbytes += write(io, UInt8('\n'), UInt8('+'))
        if has_extra_description(record)
            nbytes += unsafe_write(io, pointer(data), desclen)
        end
        # Quality
        nbytes += write(io, '\n')
        nbytes += unsafe_write(io, pointer(data) + desclen + seqlen, seqlen)
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
    print(io,   "      quality: ", collect(quality(record)))
end

function BioGenerics.isfilled(::Record)
    return true # for backwards compatibility
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end


# Accessor functions
# ------------------

function identifier(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.identifier_len)))
end

function description(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.description_len)))
end

"""
    Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)

Copy an N long block of sequence data from the fastq record `src`, starting at
position `soff`, to the `BioSequence` dest, starting at position `doff`.
"""
function Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)
    # This check is here to prevent boundserror when indexing src.sequence
    iszero(N) && return dest
    return copyto!(dest, doff, src.data, Int(src.description_len) + soff, N)
end

"""
    quality_iter(record::Record, [part::UnitRange])::Vector{UInt8}

Get an iterator of base quality of `record`. This iterator is corrupted if the record is mutated.
"""
function quality(record::Record, part::UnitRange{<:Integer}=1:seqlen(record))
    quality(record, DEFAULT_ENCODING, part)
end

function quality(record::Record, encoding::QualityEncoding, part::UnitRange{<:Integer}=1:seqlen(record))
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

Get the base quality of `record` by decoding with `encoding_name`.

The `encoding_name` can be either `:sanger`, `:solexa`, `:illumina13`, `:illumina15`, or `:illumina18`.
"""
function quality(
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
    quality(record, encoding, part)
end

function BioGenerics.seqname(record::Record)
    return identifier(record)
end

function BioGenerics.hasseqname(record::Record)
    return true
end

function BioGenerics.sequence(record::Record)
    return sequence(record)
end

function BioGenerics.sequence(::Type{S}, record::Record) where S <: BioSequences.LongSequence
    return sequence(S, record)
end

function BioGenerics.hassequence(record::Record)
    return true
end

function Base.hash(record::Record, h::UInt)
    h = hash(record.description_len, h)
    hash(view(record.data, filled(record)), h)
end
