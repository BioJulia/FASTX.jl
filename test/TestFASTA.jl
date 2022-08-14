module TestFASTA
    
using ReTest
using FASTX.FASTA: FASTA, Record, identifier, description, sequence,
    Reader, Writer, Index, index!, validate_fasta, faidx, seqlen, extract, seekrecord
using BioSequences: LongDNA, LongRNA, LongAA, @dna_str, @rna_str, @aa_str
using Random: rand!, shuffle!
using FormatSpecimens: list_valid_specimens, list_invalid_specimens, path_of_format, filename, hastag
using StringViews: StringView
using TranscodingStreams: NoopStream

@testset "Record" begin
    include("fasta/record.jl")
end
@testset "IO" begin
    include("fasta/io.jl")
end
@testset "Index" begin
    include("fasta/index.jl")
end
@testset "Specimens" begin
    include("fasta/specimens.jl")
end

end # module TestFASTA
