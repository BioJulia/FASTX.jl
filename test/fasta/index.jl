INDEX_GOOD = "abc\t100\t5\t15\t16\r\n^def?@l~2:/\t17\t1\t14\t16"

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

@testset "Parsing index" begin
    ind = Index(IOBuffer(INDEX_GOOD))
    @test ind.names["abc"] == 1
    @test ind.names["^def?@l~2:/"] == 2
    @test ind.lengths == [100, 17]
    @test ind.offsets == [5, 1]

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
    for i in 1:25
        print(buf, random_name(), '\t')
        len = rand(20:250)
        print(buf, len, '\t')
        print(buf, rand(1:100), '\t')
        lineb = rand(5:len)
        print(buf, lineb, '\t')
        linewidth = lineb + rand(1:2)
        print(buf, linewidth, '\n')
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

function make_random_indexable_fasta()
    buf = IOBuffer()
    names = String[]
    newlines = String[]
    linelengths = Int[]
    lengths = Int[]

    for i in 1:25
        newline = rand(("\n", "\r\n"))
        push!(newlines, newline)
        len = rand(20:100)
        push!(linelengths, len)
        name = random_name()
        push!(names, name)
        print(buf, '>', name, newline)
        nlines = rand(1:10)
        push!(lengths, nlines * len)
        for i in 1:nlines
            print(buf, random_seqline(len), newline)
        end
    end
    return (take!(buf), names, newlines, linelengths, lengths)
end

@testset "Creating index" begin
    @test Index(IOBuffer("")) isa Index

    for i in 1:10
        (buffer, names, newlines, linelengths, lengths) = make_random_indexable_fasta()
        index = faidx(IOBuffer(buffer))
        @test index.names == Dict(name => i for (i, name) in enumerate(names))
        @test index.lengths == lengths
        obs_lw = [FASTA.linebases_width(index, i) for i in eachindex(lengths)]
        exp_lw = [(linelengths[i], linelengths[i] + length(newlines[i])) for i in eachindex(lengths)]
        @test obs_lw == exp_lw
    end
end

@testset "Reader with index" begin
    (buffer, names, newlines, linelengths, lengths) = make_random_indexable_fasta()
    idx = faidx(IOBuffer(buffer))
    reader = Reader(IOBuffer(buffer), index=idx)

    # Test getindex
    for (i, name) in enumerate(names)
        record = reader[name]
        @test identifier(record) == name
        @test seqlen(record) == lengths[i]
    end

    # Test seekrecord
    for (i, name) in enumerate(names)
        FASTA.seekrecord(reader, name)
        record = first(reader)
        @test identifier(record) == name
        @test seqlen(record) == lengths[i]
    end

    # Test extract
    for (i, name) in enumerate(names)
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