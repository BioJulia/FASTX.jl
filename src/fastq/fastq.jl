"""
    FASTA

Module under FASTX with code related to FASTA files.
"""
module FASTQ

using Automa: Automa, @re_str, @mark, @markpos, @relpos, @abspos, onenter!, onexit!
import BioGenerics: BioGenerics
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream, NoopStream
import ..FASTX: identifier, description, sequence, seqsize, truncate, memcmp, appendfrom!, CONTEXT, throw_parser_error

const Re = Automa.RegExp

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")

end
