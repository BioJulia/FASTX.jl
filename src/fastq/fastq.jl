module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioSymbols
import BioSequences
import BioGenerics: BioGenerics, isfilled
import BioGenerics.Automa: State
import TranscodingStreams: TranscodingStreams, TranscodingStream
import ..FASTX: identifier, description, sequence

include("quality.jl")
include("record.jl")
include("readrecord.jl")
include("reader.jl")
include("writer.jl")
include("fastqread.jl")

function interleave(R1::Reader, R2::Reader, out::Writer)
    R1_buff = FASTQ.Record()
    R2_buff = FASTQ.Record()
    while !eof(R1) && !eof(R2)
        read!(x, R1_buff)
        read!(y, R2_buff)
        write(out, R1_buff)
        write(out, R2_buff)
    end
end

end
