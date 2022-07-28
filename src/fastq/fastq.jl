module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioSequences
import BioGenerics: BioGenerics
import BioGenerics.Automa: State
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..FASTX: identifier, description, sequence, UTF8, seqlen

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")

end
