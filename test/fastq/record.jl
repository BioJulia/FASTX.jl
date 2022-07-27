# Only using empty records here
@testset "Basic properties" begin
    # Equality of empty records
    record = Record()
    record2 = Record()
    @test record == record2
    @test record !== record2
    empty!(record)
    @test record == record2

    # Copying
    resize!(record.data, 1000)
    cp = copy(record)
    empty!(record.data)
    @test record == cp
    @test record.data !== cp.data

    # Components of empty records
    @test identifier(record) == ""
    @test description(record) == ""
    @test sequence(record) == ""
    @test isempty(collect(quality(record)))
end

TEST_RECORD_STRINGS = [
    # Standard records
    "@some_header\r\nAAGG\r\n+\r\njjll",
    "@prkl_19900 [a b]:211\nkjmn\n+\naabb",
    "@some_header\nAAGG\n+some_header\njjll\n\n", # same as #1

    # Edge cases:
    "@\nTAG\n+\n!!!", # empty description
    "@ ||;;211name \nkakana\n+\naabbcc", # empty some_identifier
    "@header here\n\n+\n", # empty sequence

    # 
]

TEST_BAD_RECORD_STRINGS = [
    "@some\n\nTAG\n+\r\njjj", # extra newline
    "@abc\nABC\n+\nABCD", # qual too long,
    "@abc\nABC\n+\nAB", # qual too short,
    "@A B \nC\n+A B\nA", # second header different
    "@A\nC\n+AB\nA", # second header too long
    "@AB\nC\n+A\nA", # second header too short
]

@testset "Basic construction" begin
    function test_is_equal(a::Record, b::Record)
        @test a == b
        @test identifier(a) == identifier(b)
        @test description(a) == description(b)
        @test sequence(String, a) == sequence(String, b)
        @test collect(quality(a)) == collect(quality(b))
    end

    string = "@some header\nAAGG\n+\njjll"
    record = Record(string)
    record2 = Record()

    # Identity and emptiness
    @test record != record2
    @test empty!(copy(record)) == record2

    # Test basic properties
    @test identifier(record) == "some"
    @test description(record) == "some header"
    @test sequence(record) == "AAGG"
    @test collect(quality(record)) == [Int8(i)-33 for i in "jjll"]

    # Construct from two strings and quality
    record3 = Record("some header", "AAGG", [73, 73, 75, 75])
    test_is_equal(record, record3)
    @test_throws Exception Record("some_header", "TAG", [73, 73])

    # From substrings
    record4 = Record(SubString(string, 1:lastindex(string)))
    test_is_equal(record, record4)

    # From arrays
    cu = codeunits(string)
    test_is_equal(Record(cu), record)
    test_is_equal(Record(collect(cu)), record)
end

@testset "Construction edge cases" begin
    # Can construct good examples
    strings = TEST_RECORD_STRINGS[1:6]
    records = map(strings) do string
        @test Record(string) isa Any # does not throw
        Record(string)
    end

    @test description(records[4]) == ""
    @test identifier(records[5]) == ""
    @test description(records[5]) != ""
    @test sequence(records[6]) == ""
    @test isempty(collect(quality(records[6])))

    @test records[3] == records[1]

    # Throws when constructing bad examples
    for string in TEST_BAD_RECORD_STRINGS
        @test_throws Exception Record(string)
    end
end

@testset "Equality" begin
    a, b = Record(TEST_RECORD_STRINGS[1]), Record(TEST_RECORD_STRINGS[3])
    @test a == b
    push!(a.data, 0x00)
    @test a == b

    # The descr/seq break matter
    @test Record("@AA\nT\n+\nA") != Record("@A\nAT\n+\nAT")
    # Second header does not matter
    @test Record("@AA\nT\n+\nA") == Record("@AA\nT\n+AA\nA")
end

# I.e. trailing bytes in the data field not used do not matter
@testset "Noncoding bytes" begin
    record = Record(TEST_RECORD_STRINGS[1])
    resize!(record.data, 1000)
    cp = copy(record)

    for rec in (record, cp)
        @test identifier(rec) == description(rec) == "some_header"
        @test sequence(rec) == "AAGG"
        @test collect(quality(rec)) == [73, 73, 75, 75]
    end
end

@testset "Get sequence as String/StringView" begin
    records = map(Record, TEST_RECORD_STRINGS)

    @test sequence(String, records[1]) == "AAGG"
    @test sequence(String, records[2]) == "kjmn"
    @test sequence(String, records[2], 1:3) == "kjm"
    @test sequence(String, records[2], 1:0) == ""
    @test sequence(String, records[2], 3:4) == "mn"
    @test sequence(String, records[5], 1:6) == "kakana"

    @test_throws Exception sequence(String, records[5], 0:1)
    @test_throws Exception sequence(String, records[5], 1:7)
    @test_throws Exception sequence(String, records[5], -3:3)

    # Default is StringView
    @test sequence(records[1]) isa StringView
    @test sequence(records[1]) === sequence(StringView, records[1])
    @test sequence(records[1]) == sequence(String, records[1])
    @test sequence(records[1], 2:3) == "AG"
    @test_throws Exception sequence(records[1], 0:3)
    @test_throws Exception sequence(records[1], 2:5)
end

# The same machinery as FASTA is used, and that is much more
# thoroughly tested, so I only test a few cases here
@testset "Encode BioSequences" begin
    record1 = Record("@A\nTAG\n+\nPRJ")
    record2 = Record("@A\nYJP\n+\nACG")
    
    @test sequence(LongDNA{4}, record1) == dna"TAG"
    @test sequence(LongDNA{2}, record1) == dna"TAG"
    @test sequence(LongAA, record1) == aa"TAG"
    @test sequence(LongAA, record2) == aa"YJP"

    @test_throws Exception sequence(LongRNA{2}, record1)
    @test_throws Exception sequence(LongDNA{4}, record2)
end

@testset "Hashing" begin
    records = map(Record, TEST_RECORD_STRINGS)
    @test hash(records[1]) != hash(records[2])
    @test hash(records[1]) == hash(records[3])
    @test !isequal(records[1], records[2])
    @test isequal(records[1], records[3])

    @test length(unique(records)) == length(records) - 1

    cp = copy(records[2])
    append!(cp.data, rand(UInt8, 128))
    @test hash(cp) == hash(records[2])
    @test isequal(cp, records[2])
end
