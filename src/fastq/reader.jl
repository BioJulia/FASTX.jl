# FASTQ Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    seq_transform::Union{Function, Nothing}
    record::Record
    copy::Bool
end

"""
    FASTQ.Reader(input::IO; fill_ambiguous=nothing, copy=true)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `fill_ambiguous=nothing`: fill ambiguous symbols with the given symbol
* `copy::Bool`: iterating returns fresh copies instead of the same Record
"""
function Reader(input::IO; fill_ambiguous = nothing, copy::Bool=true)
    if fill_ambiguous === nothing
        seq_transform = nothing
    else
        seq_transform = generate_fill_ambiguous(fill_ambiguous)
    end
    record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
    if !(input isa TranscodingStream)
        stream = TranscodingStreams.NoopStream(input)
        return Reader(State(stream, 1, 1, false), seq_transform, record, copy)
    else
        return Reader(State(input, 1, 1, false), seq_transform, copy)
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
    (cs, ln, f) = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum), rdr.seq_transform)
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

function generate_fill_ambiguous(symbol::BioSymbols.DNA)
    certain = map(UInt8, ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'))
    # return transform function
    return function (data, range)
        fill = convert(UInt8, convert(Char, symbol))
        for i in range
            if data[i] ∉ certain
                data[i] = fill
            end
        end
        return data
    end
end

function index!(record::Record, data::UTF8)
    stream = TranscodingStreams.NoopStream(IOBuffer(data))
    cs, linenum, found = readrecord!(stream, record, (1, 1), nothing)
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
