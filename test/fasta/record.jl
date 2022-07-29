# Only using empty records here
@testset "Basic properties" begin
    # Equality
    record = Record()
    record2 = Record()
    @test record == record2
    @test record !== record2
    empty!(record)
    @test record == record2

    # Copying
    cp = copy(record)
    @test record == cp
    @test record.data !== cp.data

    # Components of empty records
    @test identifier(record) isa AbstractString
    @test isempty(identifier(record))

    @test identifier(record) == description(record)

    @test sequence(String, record) === ""
    @test sequence(String, record, 1:0) === ""
    @test_throws BoundsError sequence(String, record, 1:1) 
end

# Parsing from strings and arrays
@testset "Basic construction" begin
    function test_is_equal(a::Record, b::Record)
        @test a == b
        @test identifier(a) == identifier(b)
        @test description(a) == description(b)
        @test sequence(String, a) == sequence(String, b)
    end

    str = ">some_identifier \tmy_description  | text\nAAT\nTA\nCCG"
    record = parse(Record, str)
    record2 = Record()

    # Identity and emptiness
    @test record != record2
    @test empty!(copy(record)) == record2

    # Basic properties
    @test identifier(record) == "some_identifier"
    @test description(record) == "some_identifier \tmy_description  | text"
    @test sequence(String, record) == "AATTACCG"

    # Construct from two strings, description and sequence
    record3 = Record("some_identifier \tmy_description  | text", "AATTACCG")
    test_is_equal(record, record3)

    # From substrings
    record4 = parse(Record, SubString(str, 1:lastindex(str)))
    test_is_equal(record, record4)

    # From arrays
    record5 = parse(Record, codeunits(str))
    test_is_equal(record, record5)
    record6 = parse(Record, collect(codeunits(str)))
    test_is_equal(record, record6)
end

@testset "Construction edge cases" begin
    # Minimal sequence
    record = parse(Record, ">\n")
    @test "" == identifier(record) == description(record) == sequence(String, record)

    # Empty identifier
    record = parse(Record, ">\tsome header\nTAGA\n\nAAG")
    @test identifier(record) == ""
    @test description(record) == "\tsome header"
    @test sequence(String, record) == "TAGAAAG"

    # Empty description
    record = parse(Record, ">\nAAG\nWpKN.\n\n")
    @test identifier(record) == description(record) == ""
    @test sequence(String, record) == "AAGWpKN."
    
    # Empty sequence
    record = parse(Record, ">ag | kop[\ta]\n\n")
    @test identifier(record) == "ag"
    @test description(record) == "ag | kop[\ta]"
    @test sequence(String, record) == ""

    # Trailing description whitespace
    record = parse(Record, ">hdr name\t \r\npkmn\naj")
    @test identifier(record) == "hdr"
    @test description(record) == "hdr name\t "
    @test sequence(String, record) == "pkmnaj"

    # Trailing sequence whitespace
    record = parse(Record, ">here\nplKn\n.\n  \t\v\n\n  \n \n")
    @test identifier(record) == description(record) == "here"
    @test sequence(String, record) == "plKn.  \t\v   "
end

@testset "Equality" begin
    record = parse(Record, ">AAG\nWpKN.\n\n")
    record2 = parse(Record, ">AAG\n\r\nWpKN.\n\n\n\r\n")
    append!(record2.data, [0x05, 0x65, 0x81])
    @test record == record2

    record3 = parse(Record, ">AA\nGWpKN.\n\n")
    @test record != record3
end

# Tests trailing bytes in data field are OK
@testset "Noncoding bytes" begin
    record = parse(Record, ">abc\nOOJM\nQQ")
    resize!(record.data, 1000)
    @test identifier(record) == description(record) == "abc"
    @test sequence(String, record) == "OOJMQQ"

    cp = copy(record)
    @test record == cp
    @test identifier(cp) == description(record) == "abc"
    @test sequence(String, cp) == "OOJMQQ"
end

# Get sequence as String/StringView
@testset "Get sequence" begin
    record = Record(codeunits("ab cAACCAAGGTTKKKMMMM"), 2, 4, 10)
    @test sequence(String, record) == "AACCAAGGTT"
    @test sequence(String, record, 1:0) == ""
    @test sequence(String, record, 1:3) == "AAC"
    @test sequence(String, record, 6:10) == "AGGTT"
    
    @test_throws Exception sequence(String, record, 6:11)
    @test_throws Exception sequence(String, record, 0:3)

    # Default: StringView
    @test sequence(record) isa StringView
    @test sequence(record) === sequence(StringView, record)
    @test sequence(record) == sequence(String, record)
    @test sequence(record, 2:6) == "ACCAA"
end

# Encode to various biosequences
@testset "Encode sequence" begin
    # Encode to LongSequence
    record = parse(Record, ">header\naAtC\nwsdNN\n\nhhH")
    @test sequence(LongDNA{4}, record) == dna"AATCWSDNNHHH"
    @test sequence(LongAA, record) == aa"AATCWSDNNHHH"
    @test_throws Exception sequence(LongDNA{2}, record)
    @test_throws Exception sequence(LongRNA{4}, record)

    # Encode empty to longsequence of any type
    record = parse(Record, ">name\n\n")
    for S in [
        LongDNA{4}, LongDNA{2}, LongRNA{4}, LongRNA{2}, LongAA
    ]
        @test sequence(S, record) == S("")
    end
end

# Includes "unique"
@testset "Hashing" begin
    records = map(i -> parse(Record, i), [
        ">A\n\n",
        ">A\nAG",
        ">AA\nG",
    ])
    # Same as previous, but with noncoding data
    push!(records, Record(codeunits("AAGGGG"), 2, 2, 1))

    @test hash(first(records)) == hash(first(records))
    @test hash(records[end]) == hash(records[end-1])
    @test isequal(records[end], records[end-1])
    @test !isequal(records[3], records[2])
    @test length(unique(records)) == length(records) - 1
end
