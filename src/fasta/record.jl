# FASTA Record
# ============

"""
    FASTA.Record

Mutable struct representing a FASTA record as parsed from a FASTA file.
The content of the record can be queried with the following functions:
`identifier`, `description`, `sequence`.

FASTA records are un-typed, i.e. they are agnostic to what kind of data they contain.

See also: [`FASTA.Reader`](@ref), [`FASTA.Writer`](@ref)

# Examples
```jldoctest
julia> rec = parse(FASTARecord, ">some header\\nTAqA\\nCC");

julia> identifier(rec)
"some"

julia> description(rec)
"some header"

julia> sequence(rec)
"TAqACC"

julia> typeof(description(rec)) == typeof(sequence(rec)) <: AbstractString
true
```
"""
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
@inline seqlen(record::Record)::Int = record.sequence_len

"""
    FASTA.Record()

Create the default FASTA record.
"""
function Record()
    return Record(Vector{UInt8}(), 0, 0, 0)
end

# This is extended with its only method in FASTX.jl, but must
# be defined here to be under the FASTA module.
function Record! end

function Base.empty!(record::Record)
    # Do not truncate the underlying data buffer
    record.identifier_len = 0
    record.description_len = 0
    record.sequence_len = 0
    return record
end

function Base.parse(::Type{Record}, data::AbstractVector{UInt8})
    # Error early on empty data to not construct buffers
    isempty(data) && throw(ArgumentError("Cannot parse empty string as FASTA record"))

    record = Record(Vector{UInt8}(undef, sizeof(data)), 0, 0, 0)
    stream = NoopStream(IOBuffer(data), bufsize=sizeof(data))
    cs, _, found = readrecord!(stream, record, (1, 1))
    
    # If found is not set, then the data terminated early
    found || throw(ArgumentError("Incomplete FASTA record"))

    # In this case, the machine ran out of data exactly after one record
    p = stream.state.buffer1.bufferpos
    p > sizeof(data) && iszero(cs) && return record

    # Else, we check all trailing data to see it contains only \r\n
    for i in p-1:sizeof(data)
        if !in(data[i], (UInt8('\r'), UInt8('\n')))
            throw(ArgumentError("Invalid trailing data after FASTA record"))
        end
    end
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
    return parse(Record, take!(buf))
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
    print(io,   "   sequence: \"", truncate(sequence(record), 40), '"')
end

# TODO: Base's hash does not hash all elements. Do we have a better implementation?
function Base.hash(record::Record, h::UInt)
    # The description length is informative of the record's content
    # in a way that the sequence length and identifier length isn't.
    # I.e. you could have ">A\nAG" vs ">AA\nG"
    h = hash(record.description_len, h)
    hash(view(record.data, filled(record)), h)
end
