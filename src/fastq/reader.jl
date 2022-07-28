# FASTQ Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    record::Record
    copy::Bool
end

"""
    FASTQ.Reader(input::IO; fill_ambiguous=nothing, copy=true)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `copy::Bool`: iterating returns fresh copies instead of the same Record
"""
function Reader(input::IO; copy::Bool=true)
    record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
    if !(input isa TranscodingStream)
        stream = TranscodingStreams.NoopStream(input)
        return Reader(State(stream, 1, 1, false), record, copy)
    else
        return Reader(State(input, 1, 1, false), record, copy)
    end
end

function Base.iterate(rdr::Reader, state=nothing)
    (cs, f) = _read!(rdr, rdr.record)
    if !f
        iszero(cs) && return nothing
        # Make sure reader's record in not invalid
        empty!(rdr.record)
        error("Unexpected error when reading FASTQ file")
    end
    return if rdr.copy
        (copy(rdr.record), nothing)
    else
        (rdr.record, nothing)
    end
end

function Base.read!(rdr::Reader, rec::Record)
    (cs, f) = _read!(rdr, rdr.record)
    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("malformed FASTQ file"))
    end    
    return rec
end

function _read!(rdr::Reader, rec::Record)
    (cs, ln, f) = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum))
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f
    return (cs, f)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.close(reader::Reader)
    if reader.state.stream isa IO
        close(reader.state.stream)
    end
    return nothing
end

function index!(record::Record, data::UTF8)
    stream = TranscodingStreams.NoopStream(IOBuffer(data))
    cs, linenum, found = readrecord!(stream, record, (1, 1))
    if !found || !allspace(stream)
        throw(ArgumentError("invalid FASTQ record"))
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
