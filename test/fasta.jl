using FASTX: FASTA
using FASTX.FASTA: Record, identifier, description, sequence, Reader, Writer

using BioSequences: LongDNA, LongRNA, LongAA, @dna_str, @rna_str, @aa_str

using Random: rand!

using FormatSpecimens: list_valid_specimens, list_invalid_specimens, path_of_format

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
    record = Record(str)
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
    record4 = Record(SubString(str, 1:lastindex(str)))
    test_is_equal(record, record4)

    # From arrays
    record5 = Record(codeunits(str))
    test_is_equal(record, record5)
    record6 = Record(collect(codeunits(str)))
    test_is_equal(record, record6)
end

@testset "Construction edge cases" begin
    # Minimal sequence
    record = Record(">\n\n")
    @test "" == identifier(record) == description(record) == sequence(String, record)

    # Empty identifier
    record = Record(">\tsome header\nTAGA\n\nAAG")
    @test identifier(record) == ""
    @test description(record) == "\tsome header"
    @test sequence(String, record) == "TAGAAAG"

    # Empty description
    record = Record(">\nAAG\nWpKN.\n\n")
    @test identifier(record) == description(record) == ""
    @test sequence(String, record) == "AAGWpKN."
    
    # Empty sequence
    record = Record(">ag | kop[\ta]\n\n")
    @test identifier(record) == "ag"
    @test description(record) == "ag | kop[\ta]"
    @test sequence(String, record) == ""

    # Trailing description whitespace
    record = Record(">hdr name\t \r\npkmn\naj")
    @test identifier(record) == "hdr"
    @test description(record) == "hdr name\t "
    @test sequence(String, record) == "pkmnaj"

    # Trailing sequence whitespace
    record = Record(">here\nplKn\n.\n  \t\v\n\n  \n \n")
    @test identifier(record) == description(record) == "here"
    @test sequence(String, record) == "plKn.  \t\v   "
end

# Tests trailing bytes in data field are OK
@testset "Noncoding bytes" begin
    record = Record(">abc\nOOJM\nQQ")
    resize!(record.data, 1000)
    @test identifier(record) == description(record) == "abc"
    @test sequence(String, record) == "OOJMQQ"

    cp = copy(record)
    @test record == cp
    @test identifier(cp) == description(record) == "abc"
    @test sequence(String, cp) == "OOJMQQ"
end

# Get sequence as String
@testset "Get sequence" begin
    record = Record(codeunits("ab cAACCAAGGTTKKKMMMM"), 2, 4, 10)
    @test sequence(String, record) == "AACCAAGGTT"
    @test sequence(String, record, 1:0) == ""
    @test sequence(String, record, 1:3) == "AAC"
    @test sequence(String, record, 6:10) == "AGGTT"
    
    @test_throws Exception sequence(String, record, 6:11)
    @test_throws Exception sequence(String, record, 0:3)
end

# Encode to various biosequences
@testset "Encode sequence" begin
    # Encode to LongSequence
    record = Record(">header\naAtC\nwsdNN\n\nhhH")
    @test sequence(LongDNA{4}, record) == dna"AATCWSDNNHHH"
    @test sequence(LongAA, record) == aa"AATCWSDNNHHH"
    @test_throws Exception sequence(LongDNA{2}, record)
    @test_throws Exception sequence(LongRNA{4}, record)

    # Encode empty to longsequence of any type
    record = Record(">name\n\n")
    for S in [
        LongDNA{4}, LongDNA{2}, LongRNA{4}, LongRNA{2}, LongAA
    ]
        @test sequence(S, record) == S("")
    end
end

@testset "Copying to LongSequence" begin
    strings = [
        "ATCGTAGTAC",         # DNA 2
        "AACGMYKATNwhdvAC",   # DNA 4
        "AUTcutUAUU",         # RNA 2
        "AUGMNmuaWUAGUC",     # RNA 4
        "AGCGGACAAC",         # DNA/RNA2
        "AHCDNnnkmaAGCNvSSW", # DNA/RNA4
        "KPLMQWDCB",          # AA
        "AKLVYhkxzX",         # AA
        "BOJarleaiilvw",      # AA
        "møjsommelig",        # Invalid
        "--m!kvLMO",          # Invalid
    ]
    seqtypes = [
        LongDNA{4},
        LongDNA{2},
        LongRNA{4},
        LongRNA{2},
        LongAA
    ]
    empty_record = Record("some content", "")
    success = false
    seq = nothing
    for seqtype in seqtypes
        short_seq = seqtype()
        long_seq = seqtype(undef, 100)
        for str in strings
            record = Record("name", str)

            # Empty sequence
            copy!(short_seq, empty_record)
            @test isempty(short_seq)
            cp = copy(long_seq)
            copyto!(long_seq, empty_record)
            @test long_seq == cp

            # Nonempty sequence
            try
                seq = seqtype(str)
                success = true
            catch error
                success = false
            end
            
            if success
                # copy! will change the size, whether smaller or larger
                empty!(short_seq)
                resize!(long_seq, 100)
                copy!(long_seq, record)
                copy!(short_seq, record)
                @test length(short_seq) == ncodeunits(str)
                @test short_seq == seq == long_seq
                
                # copyto! will error if too short...
                resize!(short_seq, ncodeunits(str) - 1)
                @test_throws BoundsError copyto!(short_seq, record)
                
                # if too long, it will leave extra symbols untouched
                resize!(long_seq, 100)
                rand!(long_seq)
                rest = long_seq[ncodeunits(str)+1:100]
                copyto!(long_seq, record)
                @test long_seq[1:ncodeunits(str)] == seq
                @test long_seq[ncodeunits(str)+1:100] == rest

                # copyto with indices.
                resize!(long_seq, ncodeunits(str))
                old = copy(long_seq)
                copyto!(long_seq, 3, record, 2, 6)
                @test old[1:2] == long_seq[1:2]
                @test long_seq[3:8] == seq[2:7]
                @test long_seq[9:end] == old[9:end]

            else
                empty!(short_seq)
                resize!(long_seq, 100)
                
                # Test both, since it must throw no matter if the exception
                # is a bounds error or an alphabet incompatibility error.
                @test_throws Exception copy!(short_seq, record)
                @test_throws Exception copyto!(short_seq, record)
                @test_throws Exception copyto!(short_seq, 1, record, 1, ncodeunits(str))
                
                @test_throws Exception copy!(long_seq, record)
                @test_throws Exception copyto!(long_seq, record)
                @test_throws Exception copyto!(long_seq, 1, record, 1, ncodeunits(str))
                
            end
        end
    end
end

# Includes "unique"
@testset "Hashing" begin
    records = map(Record, [
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

@testset "Reader basics" begin
    # Empty reader
    reader = Reader(IOBuffer(""))
    @test isnothing(iterate(reader))
    close(reader)

    # Resumable
    reader = Reader(IOBuffer(">header\nTAG\nAA\n\r\n\r\n>header2\nAAA\n\nGG\n\r\n"))
    (r, s) = iterate(reader)
    @test identifier(r) == "header"
    @test sequence(String, r) == "TAGAA"
    (r, s) = iterate(reader)
    @test identifier(r) == "header2"
    @test sequence(String, r) == "AAAGG"
    @test isnothing(iterate(reader))
    close(reader)

    # Copies on iteration
    reader = Reader(IOBuffer(">A\nG\n>A\nG"))
    records = collect(reader)
    @test first(records) !== last(records)
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
    test_writer(records, r"^>\n\n>\n\n$")

    # Does not write uncoding bytes in records
    records = [
        Record(codeunits("someheader hereAACCGGTT"), 10, 15, 3),
        Record(codeunits("fewhjlkdsjepis.."), 0, 0, 0)
    ]
    test_writer(records, r"^>someheader here\nAAC\n>\n\n")

    # Lots of records to exercise the IO a little more
    # we don't test width here, that's for later
    target_buffer = IOBuffer()
    writer_buffer = IOBuffer()
    writer = Writer(writer_buffer, 0)
    for i in 1:50
        name = join(rand('A':'z'), rand(20:30))
        name2 = join(rand('A':'z'), rand(30:40))
        descr = name * ' ' * name2
        seq = join(rand(('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'), rand(1000:2000)))
        write(target_buffer, '>', descr, '\n', seq, '\n')
        write(writer, Record(descr, seq))
    end
    flush(writer)
    writer_bytes = take!(writer_buffer)
    target_bytes = take!(target_buffer)
    close(writer)
    @test writer_bytes == target_bytes
end

@testset "Writer width" begin
    header = "some data here"
    for width in (-10, 5, 25, 50)
        for seqlen in [width-1, width, 3*width, 3*width+3, 75, 200]
            seqlen < 1 && continue
            seq = join(rand('A':'Z', seqlen))
            record = Record(header, seq)
            buf = IOBuffer()
            writer = Writer(buf, width)
            write(writer, record)
            flush(writer)
            str = String(take!(buf))
            close(writer)
            target_seq = if width < 1
                seq
            else
                join(Iterators.map(join, Iterators.partition(seq, width)), '\n')
            end
            @test str == ">" * header * '\n' * target_seq * '\n'
        end
    end
end

# Records can be written, then re-read without loss,
# except arbitrary whitespace in the sequence
@testset "Round trip" begin
    strings = [
        ">abc some 
        def
        hgi",
        ">A\n\n>A B C \nlkpo",
        "
        > here | be [dragons > 1]
        polm---GA
        --PPPLLAA
        
        >and more
        AAA",
        "",
    ]
    strings = [
        join(Iterators.map(lstrip, eachline(IOBuffer(s))), '\n')
        for s in strings
    ]
    for string in strings
        read = collect(Reader(IOBuffer(string)))
        buf = IOBuffer()
        writer = Writer(buf, 0)
        for record in read
            write(writer, record)
        end
        flush(writer)
        data = String(take!(buf))
        close(writer)
        read2 = collect(Reader(IOBuffer(string)))
        @test read == read2
    end
end

@testset "Index" begin
end

@testset "Test specimens" begin
    # All valid specimens should be read, written, re-read, and the
    # second read should be identical.
    # All invalid specimens should throw an exception
    function test_valid_specimen(path)
        try
            records = open(collect, Reader, path)
            buf = IOBuffer()
            writer = Writer(buf)
            foreach(i -> write(writer, i), records)
            flush(writer)
            data = take!(buf)
            close(writer)
            records2 = collect(Reader(IOBuffer(data)))
            issame = records == records2
            if !issame
                println("Valid format not parsed properly: $path")
            end
            @test issame
        catch e
            println("Error when parsing $path")
            @test false
        end
    end

    for specimen in list_valid_specimens("FASTA")
        path = joinpath(path_of_format("FASTA"), filename(specimen))
        test_valid_specimen(path)
    end

    for specimen in list_invalid_specimens("FASTA")
        @test_throws Exception open(collect, Reader, path)
    end
end