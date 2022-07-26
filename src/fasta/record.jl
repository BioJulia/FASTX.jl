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

"Get the indices of `data` that correspond to sequence indices `part`"
function seq_data_part(record::Record, part::AbstractUnitRange)
    start, stop = first(part), last(part)
    (start < 1 || stop > record.sequence_len) && throw(BoundsError(record, start:stop))
    Int(record.description_len) + start:Int(record.description_len) + stop
end


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
Record(data::UTF8) = parse(Record, data)

function Base.parse(::Type{Record}, data::UTF8)
    record = Record(UInt8[], 0, 0, 0)
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
    return memcmp(pointer(record1.data), pointer(record2.data), filled1) == 0
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

function BioGenerics.isfilled(::Record)
    return true
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end


# Accessor functions
# ------------------

"""
    identifier(record::Record)::StringView

Get the sequence identifier of `record`. The identifier is the header
before any whitespace. If the identifier is missing, return an empty string.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.
"""
function identifier(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.identifier_len)))
end

function BioGenerics.seqname(record::Record)
    return identifier(record)
end

function BioGenerics.hasseqname(record::Record)
    return true
end

"""
    description(record::Record)::StringView

Get the description of `record`. The description is the entire header line.
Returns an `AbstractString` view into the record.
"""
function description(record::Record)::StringView
    return StringView(view(record.data, 1:Int(record.description_len)))
end

"""
    seqview(record::Record)::AbstractVector{UInt8}

Get a view of the record as an `AbstractVector{UInt8}`.
If the record has an empty sequence, return an empty view.
"""
function seqview(record::Record)::AbstractVector{UInt8}
    view(record.data, record.description_len+1:record.sequence_len)
end

"""
    sequence_iter(T, record::Record, part)

Yields an iterator of the sequence, with elements of type `T`. `T` is constructed
through `T(Char(x))` for each byte `x`. E.g. `sequence_iter(DNA, record)`.
Mutating the record will corrupt the iterator.
"""
function sequence_iter(
    ::Type{T},
    record::Record,
    part::UnitRange{<:Integer}=(1:record.sequence_len)
) where {T <: BioSymbols.BioSymbol}
    datapart = Int(record.description_len) + first(part) : Int(record.description_len) + last(part)
    data = record.data
    checkbounds(data, datapart)
    return (T(Char(@inbounds(data[i]))) for i in datapart)
end

"""
    sequence(::Type{S}, record::Record, [part::UnitRange{Int}])::S

Get the sequence of `record`.

`S` can be either a subtype of `BioSequences.BioSequence` or `String`.
If `part` argument is given, it returns the specified part of the sequence.

!!! note
    This method makes a new sequence object every time.
    If you have a sequence already and want to fill it with the sequence
    data contained in a fasta record, you can use `Base.copyto!`.
"""
function sequence(
    ::Type{S},
    record::Record,
    part::UnitRange{Int}=1:record.sequence_len
)::S where S <: BioSequences.LongSequence
    return S(@view(record.data[seq_data_part(record, part)]))
end

function sequence(
    ::Type{String},
    record::Record,
    part::UnitRange{Int}=1:record.sequence_len
)::String
    return String(record.data[seq_data_part(record, part)])
end

function Base.copy!(dest::BioSequences.LongSequence, src::Record)
    resize!(dest, UInt(src.sequence_len))
    copyto!(dest, 1, src, 1, src.sequence_len)
end

"""
    Base.copyto!(dest::BioSequences.BioSequence, src::Record)

Copy all of the sequence data from the fasta record `src` to a biological
sequence `dest`. `dest` must have a length greater or equal to the length of
the sequence represented in the fastq record. The first n elements of `dest` are
overwritten, the other elements are left untouched.
"""
function Base.copyto!(dest::BioSequences.LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, src.sequence_len)
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

function BioGenerics.sequence(::Type{S}, record::Record) where S <: BioSequences.LongSequence
    return sequence(S, record)
end

function BioGenerics.hassequence(record::Record)
    return true
end

# TODO: Base's hash does not hash all elements. Do we have a better implementation?
function Base.hash(record::Record, h::UInt)
    # The description length is informative of the record's content
    # in a way that the sequence length and identifier length isn't.
    # I.e. you could have ">A\nAG" vs ">AA\nG"
    h = hash(record.description_len, h)
    hash(view(record.data, filled(record)), h)
end
