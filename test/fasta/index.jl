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
    @test ind.linebases == [15, 14]
    @test ind.linewidths == [16, 16]

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
function make_random_index()
    buf = IOBuffer()
    for i in 1:25
        print(buf, join(rand(VALID_INDEX_CHARS, rand(10:25))), '\t')
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

@testset "Creating index" begin
end

@testset "Reader with index" begin
    # Create dummy FASTA
    # Create index from FASTA
    # Load reader of dummy FASTA
end