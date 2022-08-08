# FASTQ Writer
# ============

"""
    FASTQ.Writer(output::IO; quality_header::Union{Nothing, Bool}=nothing)

Create a data writer of the FASTQ file format.
The writer is a `BioGenerics.IO.AbstractWriter`.
Writers take ownership of the underlying IO. Mutating or closing the underlying IO
not using the writer is undefined behaviour.
Closing the writer also closes the underlying IO.

See more examples in the FASTX documentation.

See also: [`FASTQ.Record`](@ref), [`FASTQ.Reader`](@ref)

# Arguments
* `output`: Data sink to write to
* `quality_header`: Whether to print second header on the + line. If `nothing` (default),
  check the individual `Record` objects for whether they contain a second header.

# Examples
```
julia> FASTQ.Writer(open("some_file.fq", "w")) do writer
    write(writer, record) # a FASTQ.Record
end
```
"""
struct Writer{S <: TranscodingStream} <: BioGenerics.IO.AbstractWriter
    output::S
    quality_header::UInt8 # 0x00: No, 0x01: Yes, 0x02: Same as when read
end

function Writer(io::T; quality_header::Union{Nothing, Bool}=nothing) where {T <: TranscodingStream}
    qstate = quality_header === nothing ? 0x02 : UInt8(quality_header)
    Writer{T}(io, qstate)
end

Writer(io::IO; kwargs...) = Writer(NoopStream(io); kwargs...)

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

function Base.flush(writer::Writer)
    # This is, bizarrely needed for TranscodingStreams for now.
    write(writer.output, TranscodingStreams.TOKEN_END)
    flush(writer.output)
end

function Base.write(writer::Writer, record::Record)
    output = writer.output
    n = 0
    data = record.data

    desclen = UInt(record.description_len)
    seqlength = UInt(seqlen(record))

    GC.@preserve data begin
        # Header
        n += write(output, UInt8('@'))
        n += unsafe_write(output, pointer(data), desclen)
        
        # Sequence
        n += write(output, UInt8('\n'))
        n += unsafe_write(output, pointer(data) + desclen, seqlength)

        # Second header
        n += write(output, "\n+")
        # Write description in second header if either the writer is set to do that,
        # or writer is set to look at record, and record has second header
        if writer.quality_header == 0x01 || (writer.quality_header == 0x02 && has_extra_description(record))
            n += unsafe_write(output, pointer(data), desclen)
        end

        # Quality
        n += write(output, UInt8('\n'))
        n += unsafe_write(output, pointer(data) + desclen + seqlength, seqlength)

        # Final trailing newline
        n += write(output, UInt8('\n'))
    end
    return n
end
