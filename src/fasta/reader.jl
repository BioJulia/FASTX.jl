# FASTA Reader
# ============

"""
    FASTA.Reader(input::IO; index=nothing, copy::Bool=true)

Create a buffered data reader of the FASTA file format.
The reader is a `BioGenerics.IO.AbstractReader`, a stateful iterator of `FASTA.Record`.
Readers take ownership of the underlying IO. Mutating or closing the underlying IO
not using the reader is undefined behaviour.
Closing the Reader also closes the underlying IO.

See more examples in the FASTX documentation.

See also: [`FASTA.Record`](@ref), [`FASTA.Writer`](@ref)

# Arguments
* `input`: data source
* `index`: Optional random access index (currently *fai* is supported).
  `index` can be `nothing`, a `FASTA.Index`, or an `IO` in which case an index will
  be parsed from the IO, or `AbstractString`, in which case it will be treated as a path
  to a fai file.
* `copy::Bool`: iterating returns fresh copies instead of the same Record. Set to `false`
  for improved performance, but be wary that iterating mutates records.

# Examples
```jldoctest
julia> rdr = FASTAReader(IOBuffer(">header\\nTAG\\n>another\\nAGA"));

julia> records = collect(rdr); close(rdr);

julia> foreach(println, map(identifier, records))
header
another

julia> foreach(println, map(sequence, records))
TAG
AGA
```
"""
mutable struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    stream::S
    automa_state::Int
    # set to typemin(Int) if reader uses seek, then the linenum is
    # irreversibly lost.
    encoded_linenum::Int
    index::Union{Index, Nothing}
    record::Record
    copy::Bool

    function Reader{T}(io::T, index::Union{Index, Nothing}, copy::Bool) where {T <: TranscodingStream}
        record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
        new{T}(io, 1, 1, index, record, copy)
    end
end

function Reader(io::TranscodingStream; index::Union{Index, Nothing, IO, AbstractString}=nothing, copy::Bool=true)
    index!(Reader{typeof(io)}(io, nothing, copy), index)
end

Reader(io::IO; kwargs...) = Reader(NoopStream(io); kwargs...)

"""
    index!(r::FASTA.Reader, ind::Union{Nothing, Index, IO, AbstractString})

Set the index of `r`, and return `r`.
If `ind` isa `Union{Nothing, Index}`, directly set the index to `ind`.
If `ind` isa `IO`, parse the index from the FAI-formatted IO first.
If `ind` isa `AbstractString`, treat it as the path to a FAI file to parse.

See also: [`Index`](@ref), [`FASTA.Reader`](@ref)
"""
function index! end

index!(reader::Reader, index::Union{Index, Nothing}) = (reader.index = index; reader)
index!(reader::Reader, index::Union{IO, AbstractString}) = (reader.index = Index(index); reader)

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
    (cs, f) = _read!(rdr, rec)
    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("malformed FASTA file"))
    end    
    return rec
end

function _read!(rdr::Reader, rec::Record)
    enc_linenum = rdr.encoded_linenum
    cs, ln, found = readrecord!(rdr.stream, rec, (rdr.automa_state, enc_linenum))
    rdr.automa_state = cs
    # If enc_linenum is < 0, then it was unknown when entering readrecord!,
    # and so ln is meaningless.
    enc_linenum > 0 && (rdr.encoded_linenum = ln)
    return (cs, found)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.stream
end

Base.close(reader::Reader) = close(reader.stream)

function Base.getindex(reader::Reader, name::AbstractString)
    seekrecord(reader, name)
    record = Record()
    cs, _, found = readrecord!(NoopStream(reader.stream), record, (1, 1))
    @assert cs ≥ 0 && found
    return record
end

"""
    seekrecord(reader::FASTAReader, i::Union{AbstractString, Integer})

Seek `Reader` to the `i`'th record. The next iterated record with be the `i`'th record.
`i` can be the identifier of a sequence, or the 1-based record number in the `Index`.

The `Reader` needs to be indexed for this to work.
"""
function seekrecord end

function seekrecord(reader::Reader, name::AbstractString)
    index = reader.index
    if index === nothing
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader, index.names[name])
end

function seekrecord(reader::Reader, i::Integer)
    seekrecord(reader.stream, reader.index, i)
    reader.automa_state = machine.start_state
    # Make linenum unrecoverably lost
    reader.encoded_linenum = typemin(Int)
    nothing
end

"""
    extract(reader::Reader, name::AbstractString, range::Union{Nothing, UnitRange})

Extract a subsequence given by index `range` from the sequence `named` in a
`Reader` with an index. Returns a `String`.
If `range` is nothing (the default value), return the entire sequence.
"""
function extract(
    reader::Reader,
    name::AbstractString,
    range::Union{Nothing, AbstractUnitRange{<:Integer}}=nothing
)
    # Validate it has index, and index has sequence, and range
    # is inbound
    index = reader.index
    if index === nothing
        throw(ArgumentError("no index attached"))
    end
    index_of_name = index.names[name]

    len = index.lengths[index_of_name]
    checked_range = if range !== nothing
        checkbounds(1:len, range)
        isempty(range) && return ""
        range
    else
        1:len
    end
    total_bases = length(checked_range)

    # Load all required bytes into a buffer, including newlines
    (linebases, linewidth) = linebases_width(index, index_of_name)
    len_newline = linewidth - linebases
    (start_lineind_z, start_lineoff_z) = divrem(first(checked_range) - 1, linebases)
    start_offset = start_lineind_z * linewidth + start_lineoff_z

    (stop_lineind_z, stop_lineoff_z) = divrem(last(checked_range), linebases)
    stop_offset = stop_lineind_z * linewidth + stop_lineoff_z

    until_first_newline = linebases - start_lineoff_z
    buffer = Vector{UInt8}(undef, stop_offset - start_offset)
    start_file_offset = index.offsets[index_of_name] + start_offset
    seek(reader.stream, start_file_offset)
    read!(reader.stream, buffer)

    # Now remove newlines in buffer by shifting the non-newline content
    remaining = total_bases - until_first_newline
    write_index = until_first_newline + 1
    read_index = write_index + len_newline

    while remaining > 0
        n = min(linebases, remaining)
        copyto!(buffer, write_index, buffer, read_index, n)
        write_index += n
        read_index += n + len_newline
        remaining -= n
    end
    # After having removed newlines, we shrink buffer to fit
    resize!(buffer, total_bases)

    # Now check that there are no bad bytes in our buffer
    # Note: This ByteSet must correspond to the allowed bytes in
    # the FASTA machine to ensure we can seek the same FASTA files we can read
    badpos = memchr(buffer, Val(ByteSet((UInt8('\n'), UInt8('\r'), UInt8('>')))))
    if badpos !== nothing
        error("Invalid byte in FASTA sequence line: $(buffer[badpos])")
    end

    # Return the Reader to a usable state after having messed with its
    # underlying IO, then return result
    seekrecord(reader, index_of_name)
    return String(buffer)
end
