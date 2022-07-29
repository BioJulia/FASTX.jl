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

    # The specs refer to the SAM specs, which contain this regex
    name = re"[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*"
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
        if nvector == 2
            iszero(num) && error("First offset cannot be zero in a valid FAI index")
        # Number of basepairs per line obviously cannot exceed the sequence length.
        elseif nvector == 3
            if num > vectors[1][end]
                error("Bases per line exceed sequence length on line ", string(linenum))
            end
        # Linewidth is linebases plus the length of the line terminator.
        # Since we only accept \n and \r\n as line terminator, validate
        # linewidth is linebases +1 or +2.
        elseif nvector == 4
            linebases = vectors[3][end]
            if num ∉ (linebases+1, linebases+2)
                error("Linewidth must be equal to linebases +1 or +2 at line ", string(linenum))
            end
        end
        push!(vectors[nvector], num)
        num = 0
    end,
)

ctx = Automa.CodeGenContext()
@eval function read_faidx(data::Vector{UInt8})
    start = 0
    linenum = 1
    names = Dict{String, Int}()
    num = num2 = 0
    nvector = 0
    vectors = (Int[], Int[], Int[], Int[])

    GC.@preserve data begin
        #$(Automa.generate_code(index_machine, index_actions))
        $(Automa.generate_init_code(ctx, index_machine))
        $(Automa.generate_exec_code(ctx, index_machine, index_actions))
    end

    # TODO: Rely on Automa's new error code
    if !iszero(cs)
        error("Malformed index at byte $p")
    end
    return Index(names, vectors...)
end

Index(io::IO) = read_faidx(read(io))
Index(filepath::AbstractString) = open(Index, filepath)

function Base.write(io::IO, index::Index)
    # Put names dict in a sorted array
    names = Vector{String}(undef, length(index.names))
    for (name, i) in index.names
        names[i] = name
    end
    n = 0
    for i in eachindex(names)
        n +=  write(io,
            names[i], UInt8('\t'),
            string(index.lengths[i]), UInt8('\t'),
            string(index.offsets[i]), UInt8('\t'),
            string(index.linebases[i]), UInt8('\t'),
            string(index.linewidths[i]), UInt8('\n')
        )
    end
    n
end

#=
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
=#

# TODO: Generate indexing reader

# Set the reading position of `input` to the starting position of the record `name`.
function seekrecord(input::IO, index::Index, name::AbstractString)
    n_record = index.names[name]
    offset = index.offsets[n_record]
    return seek(input, offset - 1) # compensate for > symbol
end
