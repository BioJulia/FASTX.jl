module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioSymbols
import BioSequences
import BioGenerics: BioGenerics
import BioGenerics.Automa: State
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..FASTX: identifier,
    description, header,
    sequence, UTF8, seqlen

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")
include("fastqread.jl")

end
