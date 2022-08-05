# FASTA File Format
# =================

"""
    FASTA

Module under FASTX with code related to FASTA files.
"""
module FASTA

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioGenerics: BioGenerics
import BioGenerics.Automa: State
import BioSequences
import StringViews: StringView
import TranscodingStreams: TranscodingStreams, TranscodingStream

# Trivial use, I only use it here because it's a dep of Automa anyway.
# Can be removed with no big problems
using ScanByte: memchr, ByteSet
import ..FASTX: identifier, description, sequence, UTF8, seqlen

include("record.jl")
include("readrecord.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end
