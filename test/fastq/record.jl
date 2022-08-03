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
    @test isempty(quality(record))
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
    "@AB\nC\n+A\nA", # second header too short,
    "@AB\nC\n+AB\n\t", # qual not in range
    "@AB\nABC\n+\nK V", # qual not in range
    "@AB\nABC\n+\nK\x7fV", # qual not in range
]

@testset "Basic construction" begin
    function test_is_equal(a::Record, b::Record)
        @test a == b
        @test identifier(a) == identifier(b)
        @test description(a) == description(b)
        @test sequence(String, a) == sequence(String, b)
        @test quality(a) == quality(b)
    end

    string = "@some header\nAAGG\n+\njjll"
    record = parse(Record, string)
    record2 = Record()

    # Identity and emptiness
    @test record != record2
    @test empty!(copy(record)) == record2

    # Test basic properties
    @test identifier(record) == "some"
    @test description(record) == "some header"
    @test sequence(record) == "AAGG"
    @test quality(record) == "jjll"

    # Construct from two strings and quality
    record3 = Record("some header", "AAGG", [Int8(i)-OFFSET for i in "jjll"])
    test_is_equal(record, record3)
    @test_throws Exception Record("some_header", "TAG", [Int8(i)-OFFSET for i in "jj"])

    # From substrings
    record4 = parse(Record, SubString(string, 1:lastindex(string)))
    test_is_equal(record, record4)

    # From arrays
    cu = codeunits(string)
    test_is_equal(parse(Record, cu), record)
    test_is_equal(parse(Record, collect(cu)), record)

    # From AbstractString
    record5 = parse(Record, Test.GenericString(string))
    test_is_equal(record5, record)

    # From BioSequence/quality string
    record6 = Record("some header", dna"AAGG", "jjll")
    test_is_equal(record6, record)
end

@testset "Construction edge cases" begin
    # Can construct good examples
    strings = TEST_RECORD_STRINGS[1:6]
    records = map(strings) do string
        @test parse(Record, string) isa Any # does not throw
        parse(Record, string)
    end

    @test description(records[4]) == ""
    @test identifier(records[5]) == ""
    @test description(records[5]) != ""
    @test sequence(records[6]) == ""
    @test isempty(quality(records[6]))

    @test records[3] == records[1]

    # Throws when constructing bad examples
    for string in TEST_BAD_RECORD_STRINGS
        @test_throws Exception Record(string)
    end
end

@testset "Equality" begin
    a, b = parse(Record, TEST_RECORD_STRINGS[1]), parse(Record, TEST_RECORD_STRINGS[3])
    @test a == b
    push!(a.data, 0x00)
    @test a == b

    # The descr/seq break matter
    @test parse(Record, "@AA\nT\n+\nA") != parse(Record, "@A\nAT\n+\nAT")
    # Second header does not matter
    @test parse(Record, "@AA\nT\n+\nA") == parse(Record, "@AA\nT\n+AA\nA")
end

# I.e. trailing bytes in the data field not used do not matter
@testset "Noncoding bytes" begin
    record = parse(Record, TEST_RECORD_STRINGS[1])
    resize!(record.data, 1000)
    cp = copy(record)

    for rec in (record, cp)
        @test identifier(rec) == description(rec) == "some_header"
        @test sequence(rec) == "AAGG"
        @test quality(rec) == "jjll"
    end
end

@testset "Get sequence as String/StringView" begin
    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)

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
    record1 = parse(Record, "@A\nTAG\n+\nPRJ")
    record2 = parse(Record, "@A\nYJP\n+\nACG")
    
    @test sequence(LongDNA{4}, record1) == dna"TAG"
    @test sequence(LongDNA{2}, record1) == dna"TAG"
    @test sequence(LongAA, record1) == aa"TAG"
    @test sequence(LongAA, record2) == aa"YJP"

    @test_throws Exception sequence(LongRNA{2}, record1)
    @test_throws Exception sequence(LongDNA{4}, record2)
end

# We have already tested basic quality iteration in rest of tests
@testset "Quality" begin
    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)

    # Slicing
    @test isnothing(iterate(quality(records[1], 1:0)))
    @test quality(records[1]) == "jjll"
    @test quality(records[1], 1:4) == "jjll"
    @test quality(records[1], 2:4) == "jll"

    @test_throws BoundsError quality(records[1], 0:0)
    @test_throws BoundsError quality(records[1], -2:2)
    @test_throws BoundsError quality(records[1], 4:6)

    # QualityEncoding
    @test_throws Exception QualityEncoding('B':'A', 10)
    @test_throws Exception QualityEncoding('a':'A', 10)
    @test_throws Exception QualityEncoding('Z':'Y', 10)
    @test_throws Exception QualityEncoding('A':'B', -1)
    @test_throws Exception QualityEncoding('α':'β', 10)

    # Default Quality encoding
    @test collect(quality_scores(records[2])) == [Int8(i) - OFFSET for i in "aabb"]
    @test collect(quality_scores(records[5])) == [Int8(i) - OFFSET for i in "aabbcc"]
    @test collect(quality_scores(records[2], 3:4)) == [Int8(i) - OFFSET for i in "bb"]
    @test collect(quality_scores(records[2], 4:4)) == [Int8(i) - OFFSET for i in "b"]
    
    @test_throws BoundsError quality_scores(records[2], 0:4)
    @test_throws BoundsError quality_scores(records[2], 2:5)
    @test_throws BoundsError quality_scores(records[2], 5:5)

    # Custom quality encoding
    CustomQE = QualityEncoding('A':'Z', 12)
    good = parse(Record, "@a\naaaaaa\n+\nAKPZJO")
    @test collect(quality_scores(good, CustomQE)) == [Int8(i-12) for i in "AKPZJO"]
    good = parse(Record, "@a\naaa\n+\nABC")
    @test collect(quality_scores(good, CustomQE)) == [Int8(i-12) for i in "ABC"]
    good = parse(Record, "@a\naaaa\n+\nXYZW")
    @test collect(quality_scores(good, CustomQE)) == [Int8(i-12) for i in "XYZW"]
    
    # Bad sequences
    for seq in [
        "BACDEf",
        "ABCC!",
        "abc",
        "}}!]@@"
    ]
        record = parse(Record, string("@a\n", 'a'^length(seq), "\n+\n", seq))
        @test_throws Exception collect(quality_scores(record, CustomQE))
    end

    # As strings
    @test quality(String, records[1]) == quality(StringView, records[1]) == "jjll"
    @test quality(String, records[1], 2:3) == "jl"
    @test quality(String, records[1], 1:3) == "jjl"
    @test quality(records[1]) == quality(StringView, records[1])
    @test_throws Exception quality(String, records[1], 0:3)
    @test_throws Exception quality(String, records[1], 1:5)

    @test quality(String, records[1]) isa String
    @test quality(StringView, records[1]) isa StringView
end

@testset "Named quality encodings" begin
    @test_throws Exception quality(record, :fakename)
    @test_throws Exception quality(record, :snager)
    @test_throws Exception quality(record, :illumina) # must specify which kind
    
    record = parse(Record, "@ABC\nABCDEFGHIJK\n+\n]C_Za|}~^xA")

    for (symbol, offset) in [
        (:sanger, 33),
        (:solexa, 64),
        (:illumina13, 64),
        (:illumina15, 64),
        (:illumina18, 33),
    ]
        @test collect(quality_scores(record, symbol, 1:10)) == [Int8(i - offset) for i in "]C_Za|}~^x"]
    end

    # 64-encodings do not support lower end of qual range
    bad_records = map(i -> parse(Record, i), [
        "@A\nABCDE\n+\n:KPab", # Colon not in range
        "@A\nABCDE\n+\nJKH!I", # ! not in range
        "@A\nABCDE\n+\nBDE72", # numbers not in range
    ])
    for record in bad_records
        @test_throws Exception collect(quality(record, :solexa))
        @test_throws Exception collect(quality(record, :illumina13))
        @test_throws Exception collect(quality(record, :illumina15))
    end
end

@testset "Hashing" begin
    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)
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
