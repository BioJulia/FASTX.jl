module TestFASTA

using Test    
using FASTX.FASTA: FASTA, Record, identifier, description, sequence,
    Reader, Writer, Index, index!, validate_fasta, faidx, seqsize, extract, seekrecord
using BioSequences: LongDNA, LongRNA, LongAA, @dna_str, @rna_str, @aa_str
using Random: rand!, shuffle!
using FormatSpecimens: list_valid_specimens, list_invalid_specimens, path_of_format, filename, hastag
using StringViews: StringView
using TranscodingStreams: NoopStream

const VALID_INDEX_CHARS = append!(vcat('0':'9', 'A':'Z', 'a':'z'), collect("!#\$%&+./:;?@^_|~-"))
const VALID_SEQ_BYTES = [i for i in 0x00:0xff if i ∉ UInt8.(Tuple(">\r\n"))]

include("record.jl")
include("io.jl")
include("index.jl")
include("specimens.jl")

end # module TestFASTA
