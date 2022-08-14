module TestFASTQ

# Default offset for quality
const OFFSET = 33

using ReTest
using FASTX.FASTQ: Record, Reader, Writer, identifier, description,
    sequence, quality, quality_scores, QualityEncoding, quality_header!, validate_fastq
using BioSequences: LongDNA, LongRNA, LongAA, @dna_str, @rna_str, @aa_str
using FormatSpecimens: list_valid_specimens, list_invalid_specimens, path_of_format, filename, hastag
using StringViews: StringView
using TranscodingStreams: NoopStream

@testset "Record" begin
    include("fastq/record.jl")
end
@testset "IO" begin
    include("fastq/io.jl")
end
@testset "Specimens" begin
    include("fastq/specimens.jl")
end

end # module TestFASTQ