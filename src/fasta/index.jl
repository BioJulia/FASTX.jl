# FASTA Index
# ===========
#
# Index for random access to FASTA files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    Index(src::Union{IO, AbstractString})

FASTA index object, which allows constant-time seeking of FASTA files by name.
The index is assumed to be in FAI format.

Notable methods:
* `Index(::Union{IO, AbstractString})`: Read FAI file from IO or file at path
* `write(::IO, ::Index)`: Write index in FAI format
* `faidx(::IO)::Index`: Index FASTA file
* `seekrecord(::Reader, ::AbstractString)`: Go to position of seq
* `extract(::Reader, ::AbstractString)`: Extract part of sequence

Note that the FAI specs are stricter than FASTX.jl's definition of FASTA,
such that some valid FASTA records may not be indexable.
See the specs at: http://www.htslib.org/doc/faidx.html

See also: [`FASTA.Reader`](@ref)

# Examples
```jldoctest
julia> src = IOBuffer("seqname\\t9\\t14\\t6\\t8\\nA\\t1\\t3\\t1\\t2");

julia> fna = IOBuffer(">A\\nG\\n>seqname\\nACGTAC\\r\\nTTG");

julia> rdr = FASTA.Reader(fna; index=src);

julia> seekrecord(rdr, "seqname");

julia> sequence(String, first(rdr))
"ACGTACTTG"
```
"""
struct Index
    # Vector index for the record's sequence by header: See above specification
    names::Dict{String, Int}
    lengths::Vector{Int}
    offsets::Vector{Int}
    # Upper bit is linewidth - linebases - 1, whose only valid values
    # are 0 or 1.
    encoded_linebases::Vector{UInt}

    # According to specs, the index need not be ordered by the offset in the FASTA
    # file. However, we make sure the Index object is, because it make seeking easier.
    function Index(
        names::Dict{String, Int},
        lengths::Vector{Int},
        offsets::Vector{Int},
        encoded_linebases::Vector{UInt}
    )
        issorted(offsets) && return new(names, lengths, offsets, encoded_linebases)
        perm = sortperm(offsets)
        new(
            Dict(name => perm[i] for (name, i) in names),
            lengths[perm],
            offsets[perm],
            encoded_linebases[perm]
        )
    end
end

function linebases_width(index::Index, i::Integer)
    enc = index.encoded_linebases[i]
    linebases = (enc % Int) & typemax(Int)
    linewidth = linebases + 1 + (enc ≥ (typemin(Int) % UInt))
    (linebases, linewidth)
end

function Base.show(io::IO, index::Index)
    print(io, summary(index) * ":\n")
    nrows = min(10, length(index.names))
    elided = nrows < length(index.names)
    names = Vector{String}(undef, nrows)
    found = 0
    for (name, i) in index.names
        if i ≤ nrows
            names[i] = name
            found += 1
        end
        found == nrows && break
    end
    for i in 1:nrows
        print(io, "  ")
        writeline(io, index, i, names)

        # Do not write trailing newline
        if elided || i < nrows
            write(io, UInt8('\n'))
        end
    end
    elided && print(io, '⋮')
end

index_machine = let
    newline = let
        lf = onenter!(re"\n", :countline)
        Re.opt('\r') * lf
    end

    # The specs refer to the SAM specs, which contain this regex
    name = onexit!(onenter!(re"[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*", :mark), :name)
    number = onexit!(onall!(re"[0-9]+", :digit), :number)

    line = name * re"\t" * number * re"\t" * number * re"\t" * number * re"\t" * number
    fai = Re.opt(line * Re.rep(newline * line)) * Re.rep(newline)
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
        nnum = mod1(nnum + 1, 4)
        if nnum == 1
            push!(vectors[1], num)
        elseif nnum == 2
            num < 2 && error("First offset must be at least 2 in a valid FAI index")
            push!(vectors[2], num)
        # Number of basepairs per line obviously cannot exceed the sequence length.
        elseif nnum == 3
            if num > vectors[1][end]
                error("Bases per line exceed sequence length on line ", string(linenum))
            end
            linebases = num
        # Linewidth is linebases plus the length of the line terminator.
        # Since we only accept \n and \r\n as line terminator, validate
        # linewidth is linebases +1 or +2.
        elseif nnum == 4
            if num ∉ (linebases+1, linebases+2)
                error("Linewidth must be equal to linebases +1 or +2 at line ", string(linenum))
            end
            # Encode linebases
            encoded_linebases = (linebases % UInt)
            encoded_linebases |= ifelse(num == linebases+1, UInt(0), typemin(Int) % UInt)
            push!(vectors[3], encoded_linebases)
        end
        num = 0
    end,
)

@noinline function throw_index_error(data::Vector{UInt8}, linenum::Integer, p::Integer)
    p_newline = findprev(isequal(UInt8('\n')), data, p)
    offset = p_newline === nothing ? 0 : Int(p_newline)
    col = p - offset
    error("Error when parsing FAI file: Unexpected byte at index $p (line $linenum col $col)")
end

@eval function read_faidx(data::Vector{UInt8})
    start = 0
    linenum = 1
    names = Dict{String, Int}()
    num = num2 = 0
    nnum = 0
    linebases = 0
    linebases_num = 0
    vectors = (Int[], Int[], UInt[])

    $(Automa.generate_code(CONTEXT, index_machine, index_actions))

    return Index(names, vectors...)
end

Index(io::IO) = read_faidx(read(io))
Index(filepath::AbstractString) = open(i -> Index(i), filepath)

function writeline(io::IO, index::Index, line::Integer, names::Vector{String})
    (linebases, linewidth) = linebases_width(index, line)
    write(io,
        names[line], UInt8('\t'),
        string(index.lengths[line]), UInt8('\t'),
        string(index.offsets[line]), UInt8('\t'),
        string(linebases), UInt8('\t'),
        string(linewidth),
    )
end

function Base.write(io::IO, index::Index)
    # Put names dict in a sorted array
    names = Vector{String}(undef, length(index.names))
    for (name, i) in index.names
        names[i] = name
    end
    n = 0
    for i in eachindex(names)
        n += writeline(io, index, i, names)
        n += write(io, UInt8('\n'))
    end
    n
end

function Base.print(io::IO, index::Index)
    buffer = IOBuffer()
    write(buffer, index)
    String(take!(buffer))
end

index_fasta_actions = Dict(
    :mark => :(@mark),
    :countline => :(linenum += 1),
    :identifier => quote
        identifier = unsafe_string(pointer(data, @markpos), p - @markpos)
    end,
    # Not used in fai files, the newline byte is consistent within one record
    # and since this is the first newline in a record, we set it here
    :description => quote
        uses_rn_newline = byte == UInt8('\r')
        no_more_seqlines = false

        # Disturbingly, there is no API to get the absolute position of
        # an Automa machine operating on a stream. We ought to fix this.
        # This workaround works ONLY for a NoopStream,
        # and relies on abusing the internals.
        buffer_offset = buffer.transcoded - buffer.marginpos + 1

        # We want 0-indexed, p is one-indexed, and we need the offset of first sequence
        offset = buffer_offset + p + uses_rn_newline
    end,
    :seqline => quote
        # Validate line terminator is same, i.e. no seq have have both \r\n and \n
        if p < p_end && (uses_rn_newline ⊻ (byte == UInt8('\r')))
            error("Line endings must be same within one record to index, but is not on line ", string(linenum))
        end
        # Validate sequence length is the same for all lines
        let current_seqwidth = p - @markpos
            # If on first line, seqwidth is -1, we set it correctly
            if seqwidth == -1
                seqwidth = current_seqwidth
            # If we are not supposed to see more lines, or the next line
            # is longer than expected, error
            elseif no_more_seqlines || current_seqwidth > seqwidth
                error("Sequence line width must be consistent to index, but is not on line ", string(linenum))
            # If we see a shorter line, then it must be the last line
            elseif current_seqwidth < seqwidth
                no_more_seqlines = true
            end
            seqlen += current_seqwidth
        end
    end,
    :record => quote
        record_count += 1
        
        names[identifier] = record_count
        push!(lengths, seqlen)
        push!(offsets, offset)
        enc_linebases = (seqwidth % UInt)
        enc_linebases |= ifelse(uses_rn_newline, typemin(Int) % UInt, UInt(0))
        push!(encoded_linebases, enc_linebases)

        seqwidth = -1
        seqlen = 0
    end
)

initcode = quote
    names = Dict{String, Int}()
    lengths = Int[]
    offsets = Int[]
    encoded_linebases = UInt[]

    seqwidth = -1
    seqlen = 0
    linenum = 1
    uses_rn_newline = false
    # Marks whether the current seqline must be the last in the record
    # (which is the case if its shorter than the previous)
    no_more_seqlines = false
    record_count = 0
end

returncode = quote
    if cs < 0
        throw_parser_error(data, p, linenum)
    end
    return Index(names, lengths, offsets, encoded_linebases)
end

Automa.generate_reader(
    :faidx_,
    machine,
    actions = index_fasta_actions,
    context = CONTEXT,
    initcode = initcode,
    returncode = returncode
) |> eval

# Set the reading position of `input` to the starting position of the record `name`.
function seekrecord(io::IO, index::Index, name::AbstractString)
    seekrecord(io, index, index.names[name])
end

# We seek to previous sequence to find the > start of next sequence
function seekrecord(io::IO, index::Index, i::Integer)
    i == 1 && return seekstart(io)
    linebases, linewidth = linebases_width(index, i-1)
    len = index.lengths[i-1]
    prev_offset = index.offsets[i-1]
    nlines = cld(len, linebases)
    offset = prev_offset + len + nlines * (linewidth - linebases) 
    seek(io, offset)
    return nothing
end

# Note: Current implementation relies on the IO being a NoopStream exactly,
# no other transcoding stream will do.
# This is because an indexer needs to find the absolute position of the mark
# in the stream, and this is, as far as I can tell, not supported in Automa.
# As a hacky workaround, I reach into the internals of NoopStream in the
# action dict code.

"""
    faidx(io::IO)::Index

Read a `FASTA.Index` from `io`.

See also: [`Index`](@ref)

# Examples
```jldoctest
julia> ind = faidx(IOBuffer(">ab\\nTA\\nT\\n>x y\\nGAG\\nGA"))
Index:
  ab	3	4	2	3
  x	5	14	3	4
```
"""
faidx(x::IO) = faidx_(NoopStream(x))
faidx(x::NoopStream) = faidx_(x)

# High-level interface - not sure on this yet!
"""
    faidx(fnapath::AbstractString, [idxpath::AbstractString], check=true)

Index FASTA path at `fnapath` and write index to `idxpath`.
If `idxpath` is not given, default to same name as `fnapath * ".fai"`.
If `check`, throw an error if the output file already exists

See also: [`Index`](@ref)
"""
function faidx(fnapath::AbstractString, faidxpath::AbstractString; check::Bool=true)
    check && ispath(faidxpath) && error("Output path $faidxpath already exsists")
    index = open(faidx, fnapath)
    open(i -> write(i, index), faidxpath, "w")
    index
end

faidx(path::AbstractString; check::Bool=true) = faidx(path, path * ".fai"; check=check)
