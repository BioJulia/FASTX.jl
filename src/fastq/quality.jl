# FASTQ Base Quality
# ==================
#
# A representation of positions-specific integer quality scores, as in FASTQ.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    QualityEncoding(range::StepRange{Char}, offset::Integer)

FASTQ PHRED quality encoding scheme. `QualityEncoding` objects are used to
interpret the quality scores of FASTQ records.
`range` is a range of allowed ASCII chars in the encoding, e.g. `'!':'~'` for
the most common encoding scheme.
The offset is the PHRED offset.

See also: [`quality`](@ref)

# Examples
```jldoctest
julia> read = parse(FASTQ.Record, "@hdr\nAGA\n+\nabc");

julia> qe = QualityEncoding('a':'z', 16) # hypothetical encoding

julia> collect(quality(read, qe)) == [Int8(i) - 16 for i in "abc"]
```
"""
struct QualityEncoding
    # Lowest/highest acceptable ASCII byte
    low::Int8
    high::Int8

    # ASCII byte offset, i.e. 33 for standard PHRED scores
    offset::Int8

    function QualityEncoding(ascii::StepRange{Char}, offset::Integer)
        isone(step(ascii)) || error("Must use an ordinal Char range with step 1")
        off = Int8(offset)
        (low, high) = (Int8(first(ascii)), Int8(last(ascii)))
        if low > high
            error("Quality encoding range must be nonempty")
        elseif high > 127
            error("Quality encoding only works with ASCII charsets")
        elseif offset < 0
            error("Quality offset must be non-negative")    
        elseif low < offset
            error("Low end of in quality encoding range cannot be less than offset")
        else
            return new(low, high, off)
        end
    end
end

"Sanger (Phred+33) quality score encoding"
const SANGER_QUAL_ENCODING     = QualityEncoding('!':'~', 33)

"Solexa (Solexa+64) quality score encoding"
const SOLEXA_QUAL_ENCODING     = QualityEncoding('@':'~', 64)

"Illumina 1.3 (Phred+64) quality score encoding"
const ILLUMINA13_QUAL_ENCODING = QualityEncoding('@':'~', 64)

"Illumina 1.5 (Phred+64) quality score encoding"
const ILLUMINA15_QUAL_ENCODING = QualityEncoding('B':'~', 64)

"Illumina 1.8 (Phred+33) quality score encoding"
const ILLUMINA18_QUAL_ENCODING = QualityEncoding('!':'~', 33)

const DEFAULT_ENCODING = SANGER_QUAL_ENCODING

@noinline function throw_decoding_error(encoding::QualityEncoding, ascii::Integer)
    error("Quality $ascii not in encoding range $(encoding.low):$(encoding.high)")
end

@inline function decode_quality(encoding::QualityEncoding, quality::Integer)
    check_quality(encoding, quality) || throw_decoding_error(encoding, quality)
    # We just checked it's in 0:127, so can use unsafe truncation
    (quality % Int8) - encoding.offset
end

@inline function check_quality(encoding::QualityEncoding, quality::Integer)
    (quality ≥ encoding.low) & (quality ≤ encoding.high)
end
