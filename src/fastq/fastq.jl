"""
    FASTA

Module under FASTX with code related to FASTA files.
"""
module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioSequences: BioSequences, BioSequence
import BioGenerics: BioGenerics
import BioGenerics.Automa: State
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream, NoopStream
import ..FASTX: identifier, description, sequence, UTF8, seqlen, throw_parser_error, truncate

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")

end
