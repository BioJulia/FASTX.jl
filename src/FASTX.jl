module FASTX

export
    FASTA,
    FASTQ,
    identifier,
    hasidentifier,
    description,
    hasdescription,
    header,
    sequence,
    hassequence,
    quality,
    hasquality,
    quality_iter

# Generic methods
function identifier end
function description end
function sequence end
function hasidentifier end
function hasdescription end
function hassequence end
function header end

include("fasta/fasta.jl")
include("fastq/fastq.jl")

import .FASTA
import .FASTQ
import .FASTQ: quality, hasquality, quality_iter

end # module
