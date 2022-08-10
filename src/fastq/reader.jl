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
julia> rdr = FASTQReader(IOBuffer("@readname\\nGGCC\\n+\\njk;]"));

julia> record = first(rdr); close(rdr);

julia> identifier(record)
"readname"

julia> sequence(record)
"GGCC"

julia> show(collect(quality_scores(record))) # phred 33 encoding by default
Int8[73, 74, 26, 60]
```
"""
mutable struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    stream::S
    automa_state::Int
    linenum::Int
    record::Record
    copy::Bool

    function Reader{T}(io::T, copy::Bool) where {T <: TranscodingStream}
        record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
        new{T}(io, 1, 1, record, copy)
    end
end

Reader(io::TranscodingStream; copy::Bool=true) = Reader{typeof(io)}(io, copy)
Reader(io::IO; kwargs...) = Reader(NoopStream(io); kwargs...)

function Base.iterate(rdr::Reader, state=nothing)
    (cs, f) = _read!(rdr, rdr.record)
    if !f
        iszero(cs) && return nothing
        # Make sure reader's record in not invalid
        empty!(rdr.record)
        error("Unexpected end of file when reading FASTA record")
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
    (cs, ln, found) = readrecord!(rdr.stream, rec, (rdr.automa_state, rdr.linenum))
    rdr.automa_state = cs
    rdr.linenum = ln
    return (cs, found)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

BioGenerics.IO.stream(reader::Reader) = reader.stream
Base.close(reader::Reader) = close(reader.stream)

