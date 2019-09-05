module FASTX

export
    FASTA,
    FASTQ,
    identifier,
    description,
    sequence,
    quality

# Generic methods
function identifier end
function description end
function sequence end

include("fasta/fasta.jl")
include("fastq/fastq.jl")

import .FASTA
import .FASTQ
import .FASTQ: quality

end # module
