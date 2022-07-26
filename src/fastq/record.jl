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

sequence_length(record::Record) = (record.has_description_seq_len & (typemax(Int) % UInt)) % Int
has_extra_description(record::Record) = record.has_description_seq_len ≥ (typemin(Int) % UInt)

# Number of stored bytes in data field
filled(record::Record) = record.description_len + 2 * sequence_length(record)

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
    ascii_quality = [q + offset for q in quality]
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
    seqlen = UInt(sequence_length(record))
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
    println(io)
    println(io, "  description: ", description(record))
    println(io, "     sequence: ", sequence(String, record))
    print(io,   "      quality: ", quality(record))
end

# TODO: Remove?
function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    record.quality = 1:0
    return record
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





















function Base.copy!(dest::BioSequences.LongSequence, src::Record)
    resize!(dest, seqlen(src) % UInt)
    copyto!(dest, 1, src, 1, seqlen(src))
end

"""
    Base.copyto!(dest::BioSequences.LongSequence, src::Record)

Copy all of the sequence data from the fastq record `src` to a biological
sequence `dest`. `dest` must have a length greater or equal to the length of
the sequence represented in the fastq record. The first n elements of `dest` are
overwritten, the other elements are left untouched.
"""
function Base.copyto!(dest::BioSequences.LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, seqlen(src))
end

"""
    Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)

Copy an N long block of sequence data from the fastq record `src`, starting at
position `soff`, to the `BioSequence` dest, starting at position `doff`.
"""
function Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)
    return copyto!(dest, doff, src.data, src.sequence[soff], N)
end

"""
    sequence(::Type{S}, record::Record, [part::UnitRange{Int}])

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.LongSequence` or `String`.
If `part` argument is given, it returns the specified part of the sequence.

!!! note
    This method makes a new sequence object every time.
    If you have a sequence already and want to fill it with the sequence
    data contained in a fastq record, you can use `Base.copyto!`.
"""
function sequence(::Type{S}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::S where S <: BioSequences.LongSequence
    seqpart = record.sequence[part]
    return S(@view(record.data[seqpart]))
end

"""
    sequence(::Type{String}, record::Record, [part::UnitRange{Int}])::String

Get the sequence of `record` as a String.
If `part` argument is given, it returns the specified part of the sequence.
"""
function sequence(::Type{String}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::String
    return String(record.data[record.sequence[part]])
end

"""
    sequence(record::Record, [part::UnitRange{Int}])::BioSequences.DNASequence

Get the sequence of `record` as a DNA sequence.

!!! note
    This method makes a new sequence object every time.
    If you have a sequence already and want to fill it with the sequence
    data contained in a fastq record, you can use `Base.copyto!`.
"""
function sequence(record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::BioSequences.LongDNA{4}
    return sequence(BioSequences.LongDNA{4}, record, part)
end

"Get the length of the fastq record's sequence."
@inline seqlen(record::Record) = last(record.sequence) - first(record.sequence) + 1

"""
    quality_iter(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}

Get an iterator of base quality of `record`. This iterator is corrupted if the record is mutated.
"""
function quality_iter(record::Record, offset::Integer=33, part::UnitRange{Int}=1:lastindex(record.quality))
    offs = convert(UInt8, offset)
    part = record.quality[part]
    data = record.data
    return (@inbounds(data[i]) - offs for i in part)
end

"""
    quality(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}

Get the base quality of `record`.
"""
function quality(record::Record, offset::Integer=33, part::UnitRange{Int}=1:lastindex(record.quality))::Vector{UInt8}
    collect(quality_iter(record, offset, part))
end
"""
    quality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}

Get the base quality of `record` by decoding with `encoding_name`.

The `encoding_name` can be either `:sanger`, `:solexa`, `:illumina13`, `:illumina15`, or `:illumina18`.

!!! note
    Returns `nothing` if the record has no quality string.
"""
function quality(record::Record, encoding_name::Symbol, part::UnitRange{Int}=1:lastindex(record.quality))::Vector{UInt8}
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

@deprecate hasquality(record::Record) !isempty(sequence(record))

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
	isempty(record.identifier) || (h = hash(view(record.data, record.identifier), h))
	isempty(record.description) || (h = hash(view(record.data, record.description), h))
	isempty(record.sequence) || (h = hash(view(record.data, record.sequence), h))
	isempty(record.quality) || (h = hash(view(record.data, record.quality), h))
	h
end


# Helper functions
# ----------------

function checkfilled(record)
    nothing # for backward compatibility
end
