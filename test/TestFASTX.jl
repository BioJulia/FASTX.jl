module TestFASTX

using ReTest
using FASTX: FASTA, FASTQ, identifier, description, sequence, quality
using BioSequences: LongDNA, LongAA, LongRNA
using Random: rand!

include("fastxtests.jl")
include("fasta/TestFASTA.jl")
include("fastq/TestFASTQ.jl")


end # module TestFASTX