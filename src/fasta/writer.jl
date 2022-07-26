# FASTA Writer
# ============

"""
    FASTA.Writer(output::IO; width=70)

Create a data writer of the FASTA file format.

# Arguments
* `output`: data sink
* `width=70`: wrapping width of sequence characters

# Extended help
`Writer`s take ownership of the underlying IO. That means the underlying IO may
not be directly modified such as writing or reading from it, or seeking in it.

`Writer`s carry their own buffer. This buffer is flushed when the `Writer` is closed.
Do not close the underlying IO without flushing the `Writer` first. Closing the
`Writer` automatically flushes, then closes the underlying IO, and is preferred.
"""
struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
    # maximum sequence width (no limit when width ≤ 0)
    width::Int
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

function Writer(output::IO, width::Integer=70)
    if output isa TranscodingStream
        return Writer{typeof(output)}(output, width)
    else
        stream = TranscodingStreams.NoopStream(output)
        return Writer{typeof(stream)}(stream, width)
    end
end

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
