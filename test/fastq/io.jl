@testset "Reader basics" begin
    # Empty reader
    reader = Reader(IOBuffer(""))
    @test isnothing(iterate(reader))
    close(reader)

    # Resumable
    reader = Reader(IOBuffer("@A\nTAG\n+\nJJK\n@B C\nMNB\n+B C\nLLL"))
    record = first(iterate(reader))
    @test identifier(record) == "A"
    @test sequence(record) == "TAG"
    @test collect(quality(record)) == [Int8(c - OFFSET) for c in "JJK"]
    record = first(iterate(reader))
    @test identifier(record) == "B"
    @test description(record) == "B C"
    @test sequence(record) == "MNB"
    @test collect(quality(record)) == [Int8(c - OFFSET) for c in "LLL"]
    @test isnothing(iterate(reader))
    close(reader)

    # Copies on iteration
    copy_str = "@A\nT\n+\nJ\n@A\nT\n+\nJ\n@A\nB\n+\nK"
    reader = Reader(IOBuffer(copy_str))
    records = collect(reader)
    @test records[1] == records[2]
    @test records[1] !== records[2]
    @test records[1] != records[3]
    close(reader)

    # Does not copy on iteration if copy=false
    # See comments in equivalent FASTA tests
    reader = Reader(IOBuffer(copy_str); copy=false)
    records = collect(reader)
    @test records[1] === records[2] === records[3]
    @test sequence(records[1]) == "B"
    close(reader)
end

@testset "Writer basics" begin
    function test_writer(records, regex::Regex)
        buffer = IOBuffer()
        writer = Writer(buffer)
        for record in records
            write(writer, record)
        end
        flush(writer)
        str = String(take!(buffer))
        close(writer)
        @test occursin(regex, str)
    end

    # Empty writer
    records = []
    test_writer(records, r"^$")

    # Empty records
    records = [Record(), Record()]
    test_writer(records, r"^@\n\n\+\n\n@\n\n\+\n\n$")

    # Does not write noncoding bytes in records
    records = map(Record, [
        "@ABC DEF\nkjhmn\n+ABC DEF\njjjkk",
        "@pro [1-2](HLA=2);k=1\nttagga\n+\nabcdef",
    ])
    for record in records
        resize!(record.data, 512)
    end
    test_writer(records, r"@ABC DEF\nkjhmn\n\+ABC DEF\njjjkk\n@pro \[1-2\]\(HLA=2\);k=1\nttagga\n\+\nabcdef\n")

    # Exercise the IO with larger seqs to exceed the buffer
end