module FASTX

export
    FASTA,
    FASTQ

include("fasta/fasta.jl")
include("fastq/fastq.jl")

end # module