module TestFASTQ

# Default offset for quality
const OFFSET = 33

using Test
using FASTX: FASTQ
using FASTX.FASTQ: Record, Reader, Writer, identifier, description,
    sequence, quality, quality_scores, QualityEncoding, quality_header!, validate_fastq
using BioSequences: LongDNA, LongRNA, LongAA, @dna_str, @rna_str, @aa_str
using FormatSpecimens: list_valid_specimens, list_invalid_specimens, path_of_format, filename, hastag
using StringViews: StringView
using TranscodingStreams: NoopStream

TEST_RECORD_STRINGS = [
    # Standard records
    "@some_header\r\nAAGG\r\n+\r\njjll",
    "@prkl_19900 [a b]:211\nkjmn\n+\naabb",
    "@some_header\nAAGG\n+some_header\njjll\n\n", # same as #1

    # Edge cases:
    "@\nTAG\n+\n!!!", # empty description
    "@ ||;;211name \nkakana\n+\naabbcc", # empty some_identifier
    "@header here\n\n+\n", # empty sequence

    # 
]

TEST_BAD_RECORD_STRINGS = [
    "@some\n\nTAG\n+\r\njjj", # extra newline
    "@abc\nABC\n+\nABCD", # qual too long,
    "@abc\nABC\n+\nAB", # qual too short,
    "@A B \nC\n+A B\nA", # second header different
    "@A\nC\n+AB\nA", # second header too long
    "@AB\nC\n+A\nA", # second header too short,
    "@AB\nC\n+AB\n\t", # qual not in range
    "@AB\nABC\n+\nK V", # qual not in range
    "@AB\nABC\n+\nK\x7fV", # qual not in range
]

include("record.jl")
include("io.jl")
include("specimens.jl")

end # module TestFASTQ
