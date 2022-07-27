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

    # Does not copy on iteration if copy=false
    # See comments in equivalent FASTA tests
end