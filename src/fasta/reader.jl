# FASTA Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    index::Union{Index, Nothing}
    record::Record
    copy::Bool
end

"""
    FASTA.Reader(input::IO; index=nothing, copy=true)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index`: Optional random access index (currently *fai* is supported).
  `index` can be `nothing`, a `FASTA.Index`, or an `IO` in which case an index will
  be parsed from the IO, or `AbstractString`, in which case it will be treated as a path
  to a fai file.
* `copy::Bool`: iterating returns fresh copies instead of the same Record
"""
function Reader(
    input::IO;
    index::Union{Nothing, Index, IO, AbstractString}=nothing,
    copy::Bool=true
)
    idx = if index isa Index
        index
    elseif index isa Union{AbstractString, IO}
        Index(index)
    elseif index isa Nothing
        nothing
    else
        @assert false
    end
    record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
    if !(input isa TranscodingStream)
        stream = TranscodingStreams.NoopStream(input)
        return Reader(State(stream, 1, 1, false), index, record, copy)
    else
        return Reader(State(input, 1, 1, false), index, record, copy)
    end
end

function Base.iterate(rdr::Reader, state=nothing)
    (cs, ln, f) = readrecord!(rdr.state.stream, rdr.record, (rdr.state.state, rdr.state.linenum))
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f
    if !f
        iszero(cs) && return nothing
        error("Unexpected error when reading FASTQ file")
    end
    return if rdr.copy
        (copy(rdr.record), nothing)
    else
        (rdr.record, nothing)
    end
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.read!(rdr::Reader, rec::Record)
    cs, ln, f = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum))
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f
    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("malformed FASTA file"))
    end    
    return rec
end

function Base.close(reader::Reader)
    if reader.state.stream isa IO
        close(reader.state.stream)
    end
    return nothing
end

function Base.getindex(reader::Reader, name::AbstractString)
    if reader.index === nothing
        throw(ArgumentError("no index attached"))
    end
    #seekrecord(reader.state.stream, reader.index, name)
    seekrecord(reader.state.stream, reader.index, name) # TODO: reader.index may not be type stable.
    #reader.state.cs = file_machine.start_state
    #reader.state.finished = false
    #return read(reader)
    record = Record()
    cs, linenum, found = readrecord!(TranscodingStreams.NoopStream(reader.state.stream), record, (1, 1))
    @assert cs ≥ 0 && found
    return record
end

"""
    extract(reader::Reader, Alphabet, name::AbstractString, range::UnitRange)

Extract a subsequence given by index `range` from the sequence `named` in a
`Reader` with an index. Returns a `LongSequence` with the given `Alphabet`.   
"""
function extract(reader::Reader, A::BioSequences.Alphabet, name::AbstractString, range::UnitRange)
    index = reader.index
    if index === nothing
        throw(ArgumentError("no index attached"))
    end
    i = index[name]

    # Seek to first base of range within the sequence
    offset = index.offsets[i]
    len = index.lengths[i]
    last(range) > len && throw(ArgumentError("Sequence not long enough"))
    linebase = index.linebases[i]
    linewidth = index.linewidths[i]
    len_newline = linewidth - linebase
    skip_lines = fld(first(range) - 1, linebase)
    seek(reader.state.stream, offset + (first(range) - 1) + (len_newline * skip_lines))

    # Now fill in data
    buffer = Vector{UInt8}(undef, length(range))
    filled = 0
    while filled < length(buffer)
        line = readline(reader.state.stream)
        len = min(ncodeunits(line), length(buffer) - filled)
        unsafe_copyto!(pointer(buffer, filled+1), pointer(line), len)
        filled += len
    end
    return BioSequences.LongSequence{typeof(A)}(buffer)
end

function index!(record::Record, data::UTF8)
    stream = TranscodingStreams.NoopStream(IOBuffer(data))
    cs, linenum, found = readrecord!(stream, record, (1, 1))
    if !found || !allspace(stream)
        throw(ArgumentError("invalid FASTA record"))
    end
    return record
end

function allspace(stream)
    while !eof(stream)
        if !isspace(read(stream, Char))
            return false
        end
    end
    return true
end
