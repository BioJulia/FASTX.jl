# FASTQ Reader
# ============

struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::State{S}
    seq_transform::Union{Function, Nothing}
    record::Record
    copy::Bool
end

"""
    FASTQ.Reader(input::IO; fill_ambiguous=nothing, copy=true)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `fill_ambiguous=nothing`: fill ambiguous symbols with the given symbol
* `copy::Bool`: iterating returns fresh copies instead of the same Record
"""
function Reader(input::IO; fill_ambiguous = nothing, copy::Bool=true)
    if fill_ambiguous === nothing
        seq_transform = nothing
    else
        seq_transform = generate_fill_ambiguous(fill_ambiguous)
    end
    record = Record(Vector{UInt8}(undef, 2048), 0, 0, 0)
    if !(input isa TranscodingStream)
        stream = TranscodingStreams.NoopStream(input)
        return Reader(State(stream, 1, 1, false), seq_transform, record, copy)
    else
        return Reader(State(input, 1, 1, false), seq_transform, copy)
    end
end

function Base.iterate(reader::Reader, state=nothing)
    eof(reader) && return nothing
    read!(reader, reader.record)
    return if reader.copy
        (copy(reader.record), nothing)
    else
        (reader.record, nothing)
    end
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.read!(rdr::Reader, rec::Record)
    cs, ln, f = readrecord!(rdr.state.stream, rec, (rdr.state.state, rdr.state.linenum), rdr.seq_transform)
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = f
    if !f
        cs == 0 && throw(EOFError())
        throw(ArgumentError("malformed FASTQ file"))
    end    
    return rec
end

function Base.close(reader::Reader)
    if reader.state.stream isa IO
        close(reader.state.stream)
    end
    return nothing
end

function generate_fill_ambiguous(symbol::BioSymbols.DNA)
    certain = map(UInt8, ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'))
    # return transform function
    return function (data, range)
        fill = convert(UInt8, convert(Char, symbol))
        for i in range
            if data[i] ∉ certain
                data[i] = fill
            end
        end
        return data
    end
end

function index!(record::Record, data::UTF8)
    stream = TranscodingStreams.NoopStream(IOBuffer(data))
    cs, linenum, found = readrecord!(stream, record, (1, 1), nothing)
    if !found || !allspace(stream)
        throw(ArgumentError("invalid FASTQ record"))
    end
    return record
end

function allspace(stream)
    while !eof(stream)
        if !isspace(read(stream, Char))
            return false
        end
    end
    return true
end

#=

# NOTE: This does not support line-wraps within sequence and quality.
isinteractive() && @info "Compiling FASTQ parser..."
const record_machine, file_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    rep1 = Automa.RegExp.rep1
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any
    space = Automa.RegExp.space

    hspace = re"[ \t\v]"

    header1 = let
        identifier = rep(any() \ space())
        identifier.actions[:enter] = [:mark]
        identifier.actions[:exit]  = [:header1_identifier]

        description = cat(any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:mark]
        description.actions[:exit]  = [:header1_description]

        cat('@', identifier, opt(cat(rep1(hspace), description)))
    end

    sequence = re"[A-z]*"
    sequence.actions[:enter] = [:mark]
    sequence.actions[:exit]  = [:sequence]

    header2 = let
        identifier = rep1(any() \ space())
        identifier.actions[:enter] = [:mark]
        identifier.actions[:exit]  = [:header2_identifier]

        description = cat(any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:mark]
        description.actions[:exit]  = [:header2_description]

        cat('+', opt(cat(identifier, opt(cat(rep1(hspace), description)))))
    end

    quality = re"[!-~]*"
    quality.actions[:enter] = [:mark]
    quality.actions[:exit]  = [:quality]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    record′ = cat(header1, newline, sequence, newline, header2, newline, quality)
    record′.actions[:enter] = [:anchor]
    record′.actions[:exit]  = [:record]
    record = cat(record′, newline)

    file = rep(record)

    return map(Automa.compile, (record, file))
end)()

#=
write("fastq.dot", Automa.machine2dot(file_machine))
run(`dot -Tsvg -o fastq.svg fastq.dot`)
=#

function check_identical(data1, range1, data2, range2)
    if length(range1) != length(range2) ||
       memcmp(pointer(data1, first(range1)), pointer(data2, first(range2)), length(range1)) != 0
       error("sequence and quality have non-matching header")
    end
end

function memcmp(p1::Ptr, p2::Ptr, len::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, len) % Int
end

const record_actions = Dict(
    :header1_identifier  => :(record.identifier  = (mark:p-1)),
    :header1_description => :(record.description = (mark:p-1)),
    :header2_identifier  => :(check_identical(record.data, mark:p-1, record.data, record.identifier)),
    :header2_description => :(check_identical(record.data, mark:p-1, record.data, record.description)),
    :sequence => :(record.sequence = (mark:p-1)),
    :quality  => :(record.quality  = (mark:p-1)),
    :record   => :(record.filled   = 1:p-1),
    :anchor => :(),
    :mark   => :(mark = p),
    :countline => :())
eval(
    BioCore.ReaderHelper.generate_index_function(
        Record,
        record_machine,
        :(mark = 0),
        record_actions))
eval(
    BioCore.ReaderHelper.generate_read_function(
        Reader,
        file_machine,
        :(mark = offset = 0),
        merge(record_actions, Dict(
            :header1_identifier  => :(record.identifier  = (mark:p-1) .- stream.anchor .+ 1),
            :header1_description => :(record.description = (mark:p-1) .- stream.anchor .+ 1),
            :header2_identifier  => :(check_identical(data, mark:p-1, data, (record.identifier) .+ stream.anchor .- 1)),
            :header2_description => :(check_identical(data, mark:p-1, data, (record.description) .+ stream.anchor .- 1)),
            :sequence            => :(record.sequence    = (mark:p-1) .- stream.anchor .+ 1),
            :quality             => :(record.quality     = (mark:p-1) .- stream.anchor .+ 1),
            :record => quote
                if length(record.sequence) != length(record.quality)
                    error("the length of sequence does not match the length of quality")
                end
                BioCore.ReaderHelper.resize_and_copy!(record.data, data, BioCore.ReaderHelper.upanchor!(stream):p-1)
                record.filled = (offset+1:p-1) .- offset
                if reader.seq_transform != nothing
                    reader.seq_transform(record.data, record.sequence)
                end
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor => :(BioCore.ReaderHelper.anchor!(stream, p); offset = p - 1)))))

=#
