# FASTA Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    index::Union{Index, Nothing}
end

"""
    FASTA.Reader(input::IO; index = nothing)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
function Reader(input::IO; index = nothing)
    if isa(index, AbstractString)
        index = Index(index)
    elseif index != nothing
        throw(ArgumentError("index must be a filepath or nothing"))
    end
    if !(input isa TranscodingStream)
        stream = TranscodingStreams.NoopStream(input)
    end
    return Reader(State(stream, 1, 1, false), index)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.read!(rdr::Reader, rec::Record)
    rdr.state.state == 0 && throw(EOFError())
    cs, ln, f = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum))
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f
    return rec
end

function Base.close(reader::Reader)
    if reader.state.stream isa IO
        close(reader.state.stream)
    end
    return nothing
end

function Base.getindex(reader::Reader, name::AbstractString)
    if reader.index == nothing
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

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
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
