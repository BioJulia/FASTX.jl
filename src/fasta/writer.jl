# FASTA Writer
# ============

"""
    FASTA.Writer(output::IO; width=70)

Create a data writer of the FASTA file format.
The writer is a `BioGenerics.IO.AbstractWriter`.
Writers take ownership of the underlying IO. Mutating or closing the underlying IO
not using the writer is undefined behaviour.
Closing the writer also closes the underlying IO.

See more examples in the FASTX documentation.

See also: [`FASTA.Record`](@ref), [`FASTA.Reader`](@ref)

# Arguments
* `output`: Data sink to write to
* `width`: Wrapping width of sequence characters. If < 1, no wrapping.

# Examples
```
julia> FASTA.Writer(open("some_file.fna", "w")) do writer
    write(writer, record) # a FASTA.Record
end
```
"""
mutable struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
    # maximum sequence width (no limit when width ≤ 0)
    width::Int

    function Writer{S}(output::S, width::Int) where {S <: TranscodingStream}
        finalizer(new{S}(output, width)) do writer
            close(writer.output)
        end
    end
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

Writer(io::T; width::Integer=70) where {T <: TranscodingStream} = Writer{T}(io, width)
Writer(io::IO; kwargs...) = Writer(NoopStream(io); kwargs...)

function Base.flush(writer::Writer)
    # This is, bizarrely needed for TranscodingStreams for now.
    write(writer.output, TranscodingStreams.TOKEN_END)
    flush(writer.output)
end

function Base.write(writer::Writer, record::Record)
    output = writer.output
    width = writer.width
    n::Int = 0
    # If write has no width, we can just write the record in a single line
    # as the default method does
    if width ≤ 0 || width ≥ record.sequence_len
        n += write(output, record, '\n')
    # Else we write it in chunks.
    else
        data = record.data
        GC.@preserve data begin
            # Write header
            n += write(output, UInt8('>'))
            n += unsafe_write(output, pointer(data), record.description_len)
            n += write(output, UInt8('\n'))

            # Write sequence in a loop of chunks of width bytes
            p = pointer(data, record.description_len + 1)
            p_end = pointer(data, record.description_len + record.sequence_len)
            while p ≤ p_end
                w = min(width, p_end - p + 1)
                n += unsafe_write(output, p, w)
                n += write(output, '\n')
                p += w
            end
        end
    end
    return n
end
