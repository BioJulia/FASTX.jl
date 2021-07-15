module FASTX

export
    FASTA,
    FASTQ,
    identifier,
    description,
    sequence,
    quality,
    transcode

# Generic methods
function identifier end
function description end
function sequence end

include("fasta/fasta.jl")
include("fastq/fastq.jl")

import .FASTA
import .FASTQ
import .FASTQ: quality

function FASTA.Record(record::FASTQ.Record)
    FASTQ.checkfilled(record)
    slice = (first(record.identifier) - 1):last(record.sequence)
    newdata = @inbounds record.data[slice]
    newdata[1] = 0x3E
    FASTA.Record(newdata)
end

"""
    transcode(in::FASTQ.Reader, out::FASTA.Writer)

Convert a FASTQ file to a FASTA file.
"""
function transcode(in::FASTQ.Reader, out::FASTA.Writer)
    buff_record = FASTQ.Record()
    while !eof(in)
        read!(in, buff_record)
        write(out, FASTA.Record(buff_record))
    end
end

end # module
