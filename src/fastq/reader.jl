# FASTQ Reader
# ============

"""
    FASTQ.Reader(input::IO; copy::Bool=true)

Create a buffered data reader of the FASTQ file format.
The reader is a `BioGenerics.IO.AbstractReader`, a stateful iterator of `FASTQ.Record`.
Readers take ownership of the underlying IO. Mutating or closing the underlying IO
not using the reader is undefined behaviour.
Closing the Reader also closes the underlying IO.

See more examples in the FASTX documentation.

See also: [`FASTQ.Record`](@ref), [`FASTQ.Writer`](@ref)

# Arguments
* `input`: data source
* `copy::Bool`: iterating returns fresh copies instead of the same Record. Set to `false`
  for improved performance, but be wary that iterating mutates records.

# Examples
```jldoctest
julia> rdr = Reader(IOBuffer("@readname\\nGGCC\\n+\\njk;]"));

julia> record = first(Reader); close(reader);

julia> identifier(record)
"readname"

julia> sequence(record)
"GGCC"

julia> show(collect(quality(record))) # phred 33 encoding by default
[73, 74, 26, 60]
```
"""
struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    record::Record
    copy::Bool
end

function Reader(io::TranscodingStream; copy::Bool=true)
    record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
    Reader(State(io, 1, 1, false), record, copy)
end

Reader(io::IO; kwargs...) = Reader(NoopStream(io); kwargs...)

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
    stream = NoopStream(IOBuffer(data))
    cs, linenum, found = readrecord!(stream, record, (1, 1))
    found || throw(ArgumentError("invalid FASTQ record"))
    return record
end
