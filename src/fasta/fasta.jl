# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioGenerics: BioGenerics, isfilled
import BioGenerics.Exceptions: missingerror
import BioGenerics.Automa: State
import BioSequences
import BioSymbols
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..FASTX: identifier, hasidentifier,
    description, hasdescription, header,
    sequence, hassequence

include("record.jl")
include("index.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")

end
