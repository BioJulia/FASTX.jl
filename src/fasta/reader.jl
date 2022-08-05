# FASTA Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    index::Union{Index, Nothing}
    record::Record
    copy::Bool
end

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
julia> rdr = Reader(IOBuffer(">header\nTAG\n>another\nAGA"));

julia> records = collect(Reader); close(reader);

julia> show(map(identifier, records))
["header", "another"]

julia> show(map(sequence, records))
["TAG", "AGA"]
```
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
        return Reader(State(stream, 1, 1, false), idx, record, copy)
    else
        return Reader(State(input, 1, 1, false), idx, record, copy)
    end
end

function Base.iterate(rdr::Reader, state=nothing)
    (cs, f) = _read!(rdr, rdr.record)
    if !f
        iszero(cs) && return nothing
        # Make sure reader's record in not invalid
        empty!(rdr.record)
        error("Unexpected error when reading FASTA file")
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
    cs, ln, f = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum))
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

function Base.getindex(reader::Reader, name::AbstractString)
    seekrecord(reader, name)
    record = Record()
    cs, _, found = readrecord!(TranscodingStreams.NoopStream(reader.state.stream), record, (1, 1))
    @assert cs ≥ 0 && found
    return record
end

function seekrecord(reader::Reader, name::AbstractString)
    if reader.index === nothing
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader, reader.index.names[name])
end

# TODO: `linenum` in reader.state is messed up by this operation
# It is not easily recovered.
function seekrecord(reader::Reader, i::Integer)
    seekrecord(reader.state.stream, reader.index, i)
    reader.state.state = machine.start_state
    reader.state.filled = false
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
    range::Union{Nothing, UnitRange}=nothing
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
    start_file_offset = index.offsets[index_of_name] + ncodeunits(name) + len_newline + start_offset
    seek(reader.state.stream, start_file_offset)
    read!(reader.state.stream, buffer)

    # Now remove newlines in buffer by shifting the non-newline content
    remaining = total_bases - until_first_newline
    write_index = until_first_newline + 1
    read_index = write_index + len_newline

    #=
    @show start_offset
    @show stop_offset
    @show write_index
    @show read_index
    @show len_newline
    @show remaining
    @show total_bases
    @show linebases
    =#

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

function index!(record::Record, data::UTF8)
    stream = TranscodingStreams.NoopStream(IOBuffer(data))
    _, _, found = readrecord!(stream, record, (1, 1))
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
