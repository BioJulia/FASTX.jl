INDEX_GOOD = "abc\t100\t5\t15\t16\r\n^def?@l~2:/\t17\t200\t14\t16\nABC\t12\t55\t8\t10"

INDEX_NAME_SPACE = "abc def\t100\t5\t15\t16"
INDEX_BAD_NAME_1 = "=abc\t100\t5\t15\t16"
INDEX_BAD_NAME_2 = "abc\\def\t100\t5\t15\t16"
INDEX_NEGATIVE = "abc\t100\t5\t-4\t-3"
INDEX_OVERFLOW = "abc\t10000000000000000000000000000000\t5\t15\t16"

# Longer than len
INDEX_BAD_LINEBASES = "abc\t100\t5\t101\t102"

# Too short compared to linebases
INDEX_BAD_LINEWIDTH_1 = "abc\t100\t5\t15\t15"

# Too long compared to linebases
INDEX_BAD_LINEWIDTH_2 = "abc\t100\t5\t15\t18"

INDEX_ZERO_OFFSET = "abc\t100\t5\t15\t16\ndef\t6\t0\t1\t2"

function test_same_index(a::Index, b::Index)
    @test a.names == b.names
    @test a.lengths == b.lengths
    @test a.offsets == b.offsets
    @test a.encoded_linebases == b.encoded_linebases
end

@testset "Parsing index" begin
    # Test correctly parsed _and ordered_
    ind = Index(IOBuffer(INDEX_GOOD))
    @test ind.names["abc"] == 1
    @test ind.names["ABC"] == 2
    @test ind.names["^def?@l~2:/"] == 3
    @test ind.lengths == [100, 12, 17]
    @test ind.offsets == [5, 55, 200]

    for bad_index in [
        INDEX_NAME_SPACE,
        INDEX_BAD_NAME_1,
        INDEX_BAD_NAME_2,
        INDEX_NEGATIVE,
        INDEX_OVERFLOW,
        INDEX_BAD_LINEBASES,
        INDEX_BAD_LINEWIDTH_1,
        INDEX_BAD_LINEWIDTH_2,
        INDEX_ZERO_OFFSET
    ]
        @test_throws ErrorException Index(IOBuffer(bad_index))
    end
end

const VALID_INDEX_CHARS = append!(vcat('0':'9', 'A':'Z', 'a':'z'), collect("!#\$%&+./:;?@^_|~-"))
const VALID_SEQ_BYTES = [i for i in 0x00:0xff if i ∉ UInt8.(Tuple(">\r\n"))]

random_name() = join(rand(VALID_INDEX_CHARS, rand(10:25)))
random_seqline(len::Integer) = String(rand(VALID_SEQ_BYTES, len))
function make_random_index()
    buf = IOBuffer()
    offset = 0
    for i in 1:25
        name = random_name()
        print(buf, name, '\t')
        offset += ncodeunits(name) + 1
        len = rand(20:250)
        print(buf, len, '\t')
        print(buf, offset, '\t')
        lineb = rand(5:len)
        print(buf, lineb, '\t')
        linewidth = lineb + rand(1:2)
        print(buf, linewidth, '\n')
        offset += cld(len, lineb) * (linewidth - lineb) + len
    end
    seekstart(buf)
    bytes = take!(buf)
    (Index(IOBuffer(bytes)), bytes)
end

@testset "Writing index" begin
    for i in 1:10
        index, bytes = make_random_index()
        buf = IOBuffer()
        write(buf, index)
        @test take!(buf) == bytes
    end
end

@testset "Parsing index from file" begin
    name = tempname()
    (ind, bytes) = make_random_index()
    open(name, "w") do io
        write(io, bytes)
    end
    test_same_index(ind, Index(name))
end

function make_random_indexable_fasta()
    buf = IOBuffer()
    names = String[]
    newlines = String[]
    linelengths = Int[]
    lengths = Int[]

    for i in 1:25
        newline = rand(("\n", "\r\n"))
        push!(newlines, newline)
        len = rand(100:1000)
        seq = codeunits(random_seqline(len))
        linelen = rand(10:len)
        push!(linelengths, linelen)
        name = random_name()    
        push!(names, name)
        print(buf, '>', name, newline)
        push!(lengths, len)
        for i in Iterators.partition(1:len, linelen)
            print(buf, String(seq[i]), newline)
        end
    end
    return (take!(buf), names, newlines, linelengths, lengths)
end

BADFNA_LINEENDINGS = ">abc\nTA\nAG\n>def\r\nAA\nGA\r\nGG"
BADFNA_INCONSISTENT_SEQWIDTH_1 = ">A\nTT\nTTT\nT"
BADFNA_INCONSISTENT_SEQWIDTH_2 = ">A\nTTT\nTT\nTT"

BADFNA_UNPARSEABLE = "dfklgjs\r\r\r\r\n\n\n\n"

@testset "Creating index" begin
    @test Index(IOBuffer("")) isa Index

    # Random parseable cases
    for i in 1:10
        (buffer, names, newlines, linelengths, lengths) = make_random_indexable_fasta()
        index = faidx(IOBuffer(buffer))
        @test index.names == Dict(name => i for (i, name) in enumerate(names))
        @test index.lengths == lengths
        obs_lw = [FASTA.linebases_width(index, i) for i in eachindex(lengths)]
        exp_lw = [(linelengths[i], linelengths[i] + length(newlines[i])) for i in eachindex(lengths)]
        @test obs_lw == exp_lw

        # Current implementation is special for NoopStream, test it
        index2 = faidx(NoopStream(IOBuffer(buffer)))
        test_same_index(index, index2)
    end

    # Failure cases
    for bad_case in [
        BADFNA_LINEENDINGS,
        BADFNA_INCONSISTENT_SEQWIDTH_1,
        BADFNA_INCONSISTENT_SEQWIDTH_2
    ]
        @test_throws ErrorException faidx(IOBuffer(bad_case))
    end
    @test_throws ArgumentError faidx(IOBuffer(BADFNA_UNPARSEABLE))
end

@testset "Reader with index" begin
    (buffer, names, newlines, linelengths, lengths) = make_random_indexable_fasta()
    idx = faidx(IOBuffer(buffer))
    reader1 = Reader(IOBuffer(buffer), index=idx)
    io = IOBuffer()
    write(io, idx)
    seekstart(io)
    reader2 = Reader(IOBuffer(buffer), index=io)
    name = tempname()
    open(i -> write(i, idx), name, "w")
    reader3 = Reader(IOBuffer(buffer), index=name)

    for reader in [reader1, reader2, reader3]
        # Test getindex
        inames = shuffle!(collect(enumerate(names)))

        for (i, name) in inames
            record = reader[name]
            @test identifier(record) == name
            @test seqlen(record) == lengths[i]
        end

        # Test seekrecord
        for (i, name) in inames
            FASTA.seekrecord(reader, name)
            record = first(reader)
            @test identifier(record) == name
            @test seqlen(record) == lengths[i]
        end

        # Test extract
        for (i, name) in inames
            seq = extract(reader, name)
            record = reader[name]
            @test ncodeunits(seq) == lengths[i]
            @test seq == sequence(record)

            start = rand(1:seqlen(record))
            stop = rand(start:seqlen(record))
            seq = extract(reader, name, start:stop)
            seq2 = sequence(record, start:stop)
            @test seq == seq2
        end
    end

    # Test extract with bad seq
    data = Vector{UInt8}(">A\nABC")
    reader = Reader(IOBuffer(data), index=faidx(IOBuffer(data)))
    data[5] = UInt8('>')
    @test_throws Exception extract(reader, "A")

    # Test can't do these operations without an index
    data = ">A\nT"
    reader = Reader(IOBuffer(data))
    @test_throws Exception reader["A"]
    @test_throws Exception extract(reader, "A")
end

@testset "Faidx existing file" begin
    name1 = tempname()
    name2 = tempname()
    (buffer, names, newlines, linelengths, lengths) = make_random_indexable_fasta()
    open(i -> write(i, buffer), name1, "w")

    # Generic
    faidx(name1)
    @test Index(name1 * ".fai") isa Index

    # Already exists
    @test_throws Exception faidx(name1)

    faidx(name1, name2)
    @test Index(name2) isa Index

    # Force overwrite file
    open(i -> print(i, "bad data"), name2, "w")
    @test_throws Exception faidx(name1, name2)
    faidx(name1, name2, check=false)
    @test Index(name2) isa Index
end

@testset "Issue" begin
    data = """>seq1 sequence
TAGAAAGCAA
TTAAAC
>seq2 sequence
AACGG
UUGC
"""
    reader = FASTA.Reader(IOBuffer(data), index=faidx(IOBuffer(data)))
    seekrecord(reader, "seq2")
    record = first(reader)
    @test sequence(record) == "AACGGUUGC"

    record = reader["seq1"]
    @test description(record) == "seq1 sequence"

    @test extract(reader, "seq1", 3:5) == "GAA"
end