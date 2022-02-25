# FASTQ Writer
# ============

struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
    quality_header::Bool
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

Writer(output::IO; quality_header::Bool=false) = Writer(output, quality_header)
"""
    FASTQ.Writer(output::IO; quality_header=false)

Create a data writer of the FASTQ file format.

# Arguments
* `output`: data sink
* `quality_header=false`: output the title line at the third line just after '+'

# Extended help
`Writer`s take ownership of the underlying IO. That means the underlying IO may
not be directly modified such as writing or reading from it, or seeking in it.

`Writer`s carry their own buffer. This buffer is flushed when the `Writer` is closed.
Do not close the underlying IO without flushing the `Writer` first. Closing the
`Writer` automatically flushes, then closes the underlying IO, and is preferred.
"""
function Writer(output::IO, quality_header::Bool)
    if output isa TranscodingStream
        return Writer{typeof(output)}(output, quality_header)
    else
        stream = TranscodingStreams.NoopStream(output)
        return Writer{typeof(stream)}(stream, quality_header)
    end
end

function Base.write(writer::Writer, record::Record)
    checkfilled(record)
    output = writer.output
    n = 0
    # sequence
    n += write(output, '@')
    n += unsafe_write(output, pointer(record.data, first(record.identifier)), length(record.identifier))
    if hasdescription(record)
        n += write(output, ' ')
        n += unsafe_write(output, pointer(record.data, first(record.description)), length(record.description))
    end
    n += write(output, '\n')
    n += unsafe_write(output, pointer(record.data, first(record.sequence)), length(record.sequence))
    n += write(output, '\n')
    # quality
    n += write(output, '+')
    if writer.quality_header
        n += unsafe_write(output, pointer(record.data, first(record.identifier)), length(record.identifier))
        if hasdescription(record)
            n += write(output, ' ')
            n += unsafe_write(output, pointer(record.data, first(record.description)), length(record.description))
        end
    end
    n += write(output, '\n')
    n += unsafe_write(output, pointer(record.data, first(record.quality)), length(record.quality))
    n += write(output, '\n')
    return n
end
