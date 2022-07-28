# FASTA Index
# ===========
#
# Index for random access to FASTA files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# http://www.htslib.org/doc/faidx.html
"""
    Index(src::Union{IO, AbstractString})

FASTA index object, which allows constant-time seeking of FASTA files by name.
The index is assumed to be in FAI format.

Notable methods:
* `Index(::Union{IO, AbstractString})`: Read FAI file from IO or file at path
* `write(::IO, ::Index)`: Write index in FAI format
* `faidx(::IO)::Index`: Index FASTA file

See also: [FASTA.Reader](@ref)

# Examples
```jldoctest
julia> src = IOBuffer("seqname\t9\t0\t6\t8");

julia> fna = IOBuffer(">A\nG\n>seqname\nACGTAC\r\nTTG");

julia> rdr = FASTA.Reader(fna; index=src)

julia> seekrecord(rdr, "seqname");

julia> sequence(String, first(rdr))
"ACGTACTTG"
```
"""
struct Index
    # offset for the record's sequence by header: See above specification
    names::Dict{String, Int}
    lengths::Vector{Int}
    offsets::Vector{Int}
    linebases::Vector{Int}
    linewidths::Vector{Int}
end

index_machine = let
    re = Automa.RegExp
    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        re.opt('\r') * lf
    end

    name = re.rep1(re.any() \ re.space())
    name.actions[:enter] = [:mark]
    name.actions[:exit] = [:name]

    number = re"[1-9][0-9]*"
    number.actions[:all] = [:digit]
    number.actions[:exit] = [:number]

    line = name * re"\t" * number * re"\t" * number * re"\t" * number * re"\t" * number
    fai = re.opt(line) * re.rep(newline * line) * re.rep(newline)
    Automa.compile(fai)
end

index_actions = Dict{Symbol, Expr}(
    :mark => :(start = p),
    :countline => :(linenum += 1),
    :name => quote
        let n = p - start
            name = unsafe_string(pointer(data, start), n)
            names[name] = linenum
        end
    end,
    :digit => quote
        num2 = 10*num + (byte - 0x30)
        if num2 < num
            error("Integer overflow on line " * string(linenum))
        end
        num = num2
    end,
    :number => quote
        nvector = mod1(nvector + 1, 4)
        push!(@inbounds (vectors[nvector]), num)
        num = 0
    end,
)

@eval function read_faidx(data::Vector{UInt8})
    start = 0
    linenum = 1
    names = Dict{String, Int}()
    num = num2 = 0
    nvector = 0
    vectors = (Int[], Int[], Int[], Int[])

    GC.@preserve data begin
        $(Automa.generate_code(index_machine, index_actions))
    end

    return Index(names, vectors...)
end

Index(io::IO) = read_faidx(read(io))
Index(filepath::AbstractString) = open(Index, filepath)

index_fasta_actions = Dict(
    :mark => :(@mark),
    :countline => :(linenum += 1),
    :identifier => quote
        # markpos start at first byte after >, one-indexed, we want 0-indexed offset
        identifier_offset = @abspos(@markpos - 1)
        identifier = unsafe_string(data, @markpos, @relpos(p))
    end,
    # Not used in fai files
    :description => quote nothing end,
    :seqline => quote
        # Validate line terminator is same, i.e. no seq have have both \r\n and \n
        if newline_byte == 0xff
            newline_byte = byte
        elseif newline_byte != byte
            error("Line endings must be same within one record to index, but is not on line ", string(linenum))
        end
        # Validate sequence length is the same for all lines
        let current_seqwith = @relpos(p-1)
            if seqwidth == -1
                seqwidth = current_seqwith
            elseif seqwidth != current_seqwith
                error("Sequence line width must be consistent to index, but is not on line ", string(linenum))
            end
        end
        seqlen += seqwidth
    end,
    :record => quote
        names[identifier] = record_count
        push!(lengths, seqlen)
        push!(offsets, identifier_offset)
        push!(linebases, seqwidth)
        linewidth = seqwidth + 1 + (newline_byte == UInt8('\r'))
        push!(linewidths, linewidth)

        linewidth = -1
        newline_byte = 0xff
        record_count += 1
        length = 0

        @escape
    end
)

# TODO: Generate indexing reader

# Set the reading position of `input` to the starting position of the record `name`.
function seekrecord(input::IO, index::Index, name::AbstractString)
    i = index[name]

    # For the first index, offset is trivially zero
    if i == 1
        offset = 0
    # Else, we go to the previous one and calculate the length of the previous
    # sequence in bytes, then seek to right after that one.
    else
        prev_offset = index.offsets[i - 1]
        prev_len = index.lengths[i - 1]
        prev_linebase = index.linebases[i - 1]
        prev_linewidth = index.linewidths[i - 1]

        # Note: newline_len may differ between sequences in the same file, as per
        # the specification.
        newline_len = prev_linewidth - prev_linebase
        len = cld(prev_len, prev_linebase) * newline_len + prev_len
        offset = prev_offset + len
    end
    seek(input, offset)
    return offset
end
