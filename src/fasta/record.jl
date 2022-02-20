# FASTA Record
# ============

mutable struct Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # indexes
    identifier::UnitRange{Int}
    description::UnitRange{Int}
    sequence::UnitRange{Int}
end

"""
    FASTA.Record()

Create an unfilled FASTA record.
"""
function Record()
    return Record(UInt8[], 1:0, 1:0, 1:0, 1:0)
end

"""
    FASTA.Record(data::Vector{UInt8})

Create a FASTA record object from `data`.

This function verifies and indexes fields for accessors.

!!! warning
    Note that the ownership of `data` is transferred to a new record object.
    Editing the input data will edit the record, and is not advised after
    construction of the record.
"""
function Record(data::Vector{UInt8})
    record = Record(data, 1:0, 1:0, 1:0, 1:0)
    index!(record)
    return record
end

"""
    FASTA.Record(str::AbstractString)

Create a FASTA record object from `str`.

This function verifies and indexes fields for accessors.
"""
Record(str::AbstractString) = Record(Vector{UInt8}(str))

Base.parse(::Record, str::AbstractString) = Record(str)

"""
    FASTA.Record(identifier, sequence)

Create a FASTA record object from `identifier` and `sequence`.
"""
function Record(identifier::AbstractString, sequence)
    return Record(identifier, nothing, sequence)
end

"""
    FASTA.Record(identifier, description, sequence)

Create a FASTA record object from `identifier`, `description` and `sequence`.
"""
function Record(identifier::AbstractString, description::Union{AbstractString,Nothing}, sequence::Union{BioSequences.BioSequence, AbstractString})
    buf = IOBuffer()
    print(buf, '>', strip(identifier))
    if description != nothing
        print(buf, ' ', description)
    end
    print(buf, '\n')
    print(buf, sequence, '\n')
    return Record(take!(buf))
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.identifier,
        record.description,
        record.sequence)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   identifier: ", hasidentifier(record) ? identifier(record) : "<missing>")
        println(io, "  description: ", hasdescription(record) ? description(record) : "<missing>")
          print(io, "     sequence: ", hassequence(record) ? truncate(sequence(String, record), 40) : "<missing>")
    else
        print(io, " <not filled>")
    end
end

function truncate(s::String, len::Integer)
    if length(s) > len
        return "$(String(collect(Iterators.take(s, len - 1))))…"
    else
        return s
    end
end

function initialize!(record::Record)
    record.filled = 1:0
    record.identifier = 1:0
    record.description = 1:0
    record.sequence = 1:0
    return record
end

function BioGenerics.isfilled(record::Record)
    return !isempty(record.filled)
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
    return StringView(view(record.data, record.identifier))
end

"""
    hasidentifier(record::Record)

Checks whether or not the `record` has an identifier.
"""
function hasidentifier(record)
    return isfilled(record) && !isempty(record.identifier)
end

function BioGenerics.seqname(record::Record)
    return identifier(record)
end

function BioGenerics.hasseqname(record::Record)
    return hasidentifier(record)
end

"""
    description(record::Record)::Union{StringView, Nothing}

Get the description of `record`.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.

!!! note
    Returns `nothing` if record has no description.
"""
function description(record::Record)::Union{String, Nothing}
    checkfilled(record)
    if !hasdescription(record)
        return nothing
    end
    return StringView(view(record.data, record.description))
end

"""
    hasdescription(record::Record)

Checks whether or not the `record` has a description.
"""
function hasdescription(record::Record)
    return isfilled(record) && !isempty(record.description)
end

"""
    header(record::Record)::StringView

Returns the stripped header line of `record`. If the header is empty, return an empty string.
Returns an `AbstractString` view into the record. If the record is overwritten,
the string data will be corrupted.
"""
function header(record::Record)::StringView
    id, de = record.identifier, record.description
    range = if isempty(id) && isempty(de)
        1:0
    elseif isempty(de)
        id
    elseif isempty(id)
        de
    else
        first(id):last(de)
    end
    return StringView(view(record.data, range))
end

"""
    sequence_iter(T, record::Record)

Yields an iterator of the sequence, with elements of type `T`. `T` is constructed
through `T(Char(x))` for each byte `x`. E.g. `sequence_iter(DNA, record)`.
Mutating the record will corrupt the iterator.
"""
function sequence_iter(::Type{T}, record::Record,
    part::UnitRange{<:Integer}=1:lastindex(record.sequence)) where {T <: BioSymbols.BioSymbol}
    checkfilled(record)
    seqpart = record.sequence[part]
    data = record.data
    return (T(Char(@inbounds (data[i]))) for i in seqpart)
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
function sequence(::Type{S}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::S where S <: BioSequences.LongSequence
    checkfilled(record)
    seqpart = record.sequence[part]
    return S(@view(record.data[seqpart]))
end

function sequence(::Type{String}, record::Record, part::UnitRange{Int}=1:lastindex(record.sequence))::String
    checkfilled(record)
    return String(record.data[record.sequence[part]])
end

"""
    hassequence(record::Record)

Checks whether or not a sequence record contains a sequence.
"""
function hassequence(record::Record)
    # zero-length sequence may exist
    return isfilled(record)
end

"Get the length of the fasta record's sequence."
@inline seqlen(record::Record) = last(record.sequence) - first(record.sequence) + 1

function Base.copy!(dest::BioSequences.LongSequence, src::Record)
    resize!(dest, seqlen(src) % UInt)
    copyto!(dest, 1, src, 1, seqlen(src))
end

"""
    Base.copyto!(dest::BioSequences.BioSequence, src::Record)

Copy all of the sequence data from the fasta record `src` to a biological
sequence `dest`. `dest` must have a length greater or equal to the length of
the sequence represented in the fastq record. The first n elements of `dest` are
overwritten, the other elements are left untouched.
"""
function Base.copyto!(dest::BioSequences.LongSequence, src::Record)
    return copyto!(dest, 1, src, 1, seqlen(src))
end

"""
    Base.copyto!(dest::BioSequences.BioSequence, doff, src::Record, soff, N)

Copy an N long block of sequence data from the fasta record `src`, starting at
position `soff`, to the `BioSequence` dest, starting at position `doff`.
"""
function Base.copyto!(dest::BioSequences.LongSequence, doff, src::Record, soff, N)
    checkfilled(src)
    if !hassequence(src)
        missingerror(:sequence)
    end
    
    # This check is here to prevent boundserror when indexing src.sequence
    iszero(N) && return dest
    return copyto!(dest, doff, src.data, src.sequence[soff], N)
end

function BioGenerics.sequence(::Type{S}, record::Record) where S <: BioSequences.LongSequence
    return sequence(S, record)
end

function BioGenerics.hassequence(record::Record)
    return hassequence(record)
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled FASTA record"))
    end
end

# Predict sequence type based on character frequencies in `seq[start:stop]`.
function predict_seqtype(seq::Vector{UInt8}, range)
    # count characters
    a = c = g = t = u = n = alpha = 0
    for i in range
        @inbounds x = seq[i]
        if x == 0x41 || x == 0x61
            a += 1
        elseif x == 0x43 || x == 0x63
            c += 1
        elseif x == 0x47 || x == 0x67
            g += 1
        elseif x == 0x54 || x == 0x74
            t += 1
        elseif x == 0x55 || x == 0x75
            u += 1
        elseif x == 0x4e || x == 0x6e
            n += 1
        end
        if 0x41 ≤ x ≤ 0x5a || 0x61 ≤ x ≤ 0x7a
            alpha += 1
            if alpha ≥ 300 && t + u > 0 && a + c + g + t + u + n == alpha
                # pretty sure that the sequence is either DNA or RNA
                break
            end
        end
    end

    # the threshold (= 0.95) is somewhat arbitrary
    if (a + c + g + t + u + n) / alpha > 0.95
        if t ≥ u
            return BioSequences.LongDNA{4}
        else
            return BioSequences.LongRNA{4}
        end
    else
        return BioSequences.LongAA
    end
end

function Base.hash(record::Record, h::UInt)
	isempty(record.identifier) || (h = hash(view(record.data, record.identifier), h))
	isempty(record.description) || (h = hash(view(record.data, record.description), h))
	isempty(record.sequence) || (h = hash(view(record.data, record.sequence), h))
	h
end
