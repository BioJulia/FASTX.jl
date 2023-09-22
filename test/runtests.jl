module FASTXTests

export TestFASTA
export TestFASTQ

using FASTX: FASTA, FASTQ, identifier, description, sequence, quality
using BioSequences: LongDNA, LongAA, LongRNA
using Random: rand!
using Test

include("maintests.jl")
include("fasta/TestFASTA.jl")
include("fastq/TestFASTQ.jl")

end # module FASTXTests
