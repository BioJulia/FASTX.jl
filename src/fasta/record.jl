# FASTA Record
# ============

mutable struct Record
    # Data contains the description, then the sequence immediately after
    # without newlines, or the initial > symbol, and then any unused trailing bytes
    data::Vector{UInt8}

    # Identifier is data[1:identifier_len]
    identifier_len::Int32

    # Description is data[1:description_len], i.e. it includes the identifier
    description_len::Int32

    # Sequence is data[description_len+1 : description_len+sequence_len]
    sequence_len::Int
end

filled(x::Record) = Int(x.description_len) + Int(x.sequence_len)
@inline seqlen(record::Record) = record.sequence_len

"""
    FASTA.Record()

Create the default FASTA record.
"""
function Record()
    return Record(Vector{UInt8}(), 0, 0, 0)
end

function Base.empty!(record::Record)
    # Do not truncate the underlying data buffer
    record.identifier_len = 0
    record.description_len = 0
    record.sequence_len = 0
    return record
end

"""
    FASTA.Record(data::Union{Vector{UInt8, String, SubString{String}}})

Create a FASTA record object from `data`.

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

"""
    FASTA.Record(description::AbstractString, sequence)

Create a FASTA record object from `description` and `sequence`.
"""
function Record(description::AbstractString, sequence::Union{BioSequences.BioSequence, AbstractString})
    buf = IOBuffer()
    print(buf, '>', description, '\n')
    # If the sequence is empty, we need to print a newline in order to not
    # have the FASTA file truncated, thus invalid
    print(buf, isempty(sequence) ? '\n' : sequence)
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
        record.sequence_len
    )
end

function Base.write(io::IO, record::Record)
    data = record.data
    write(io, UInt8('>'))
    GC.@preserve data begin
        unsafe_write(io, pointer(data), UInt(record.description_len))
        write(io, UInt8('\n'))
        unsafe_write(io, pointer(data) + UInt(record.description_len), UInt(record.sequence_len))
    end
    return filled(record) + 2 # number of bytes
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    println(io)
    println(io, "description: \"", description(record), '"')
    print(io,   "   sequence: \"", truncate(sequence(String, record), 40), '"')
end

function truncate(s::String, len::Integer)
    if length(s) > len
        return string(String(collect(Iterators.take(s, len - 1))), '…')
    else
        return s
    end
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
    Base.copyto!(dest::BioSequences.BioSequence, doff, src::Record, soff, N)

Copy an N long block of sequence data from the fasta record `src`, starting at
position `soff`, to the `BioSequence` dest, starting at position `doff`.
"""
function Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)
    # This check is here to prevent boundserror when indexing src.sequence
    iszero(N) && return dest
    return copyto!(dest, doff, src.data, Int(src.description_len) + soff, N)
end

# TODO: Base's hash does not hash all elements. Do we have a better implementation?
function Base.hash(record::Record, h::UInt)
    # The description length is informative of the record's content
    # in a way that the sequence length and identifier length isn't.
    # I.e. you could have ">A\nAG" vs ">AA\nG"
    h = hash(record.description_len, h)
    hash(view(record.data, filled(record)), h)
end
