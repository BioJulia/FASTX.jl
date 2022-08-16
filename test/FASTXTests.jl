module FASTXTests

export TestFASTA
export TestFASTQ

using ReTest
using FASTX: FASTA, FASTQ, identifier, description, sequence, quality
using BioSequences: LongDNA, LongAA, LongRNA
using Random: rand!

include("maintests.jl")
include("fasta/TestFASTA.jl")
include("fastq/TestFASTQ.jl")


end # module FASTXTests