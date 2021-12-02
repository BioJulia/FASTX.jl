module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioSymbols
import BioSequences
import BioGenerics: BioGenerics, isfilled
import BioGenerics.Automa: State
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..FASTX: identifier, hasidentifier,
    description, hasdescription, header,
    sequence, hassequence

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")
include("fastqread.jl")

end
