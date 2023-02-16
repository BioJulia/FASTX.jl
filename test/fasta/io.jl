@testset "IO" begin


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
    copy_str = ">A\nG\n>A\nG\n>A\nT"
    reader = Reader(IOBuffer(copy_str))
    records = collect(reader)
    @test records[1] == records[2]
    @test records[1] !== records[2]
    @test records[1] != records[3]
    close(reader)

    # Test Base.read! works
    reader = Reader(IOBuffer(">header string\r\nYWBL\nKKL\r\n>another\nAAGTC"))
    record = Record()
    read!(reader, record)
    @test identifier(record) == "header"
    @test description(record) == "header string"
    @test sequence(record) == "YWBLKKL"
    (record, _) = iterate(reader)
    @test (description(record), sequence(record)) == ("another", "AAGTC")

    # Does not copy on iteration if copy=false
    # in this case it will iterate the same object which
    # will just be overwritten.
    # This is not intended to be relied on, so can be removed,
    # but is currently necessary for performance, so we test it
    # to prevent performance regressions on this front
    reader = Reader(IOBuffer(copy_str); copy=false)
    records = [first(iterate(reader)) for i in 1:3]
    @test records[1] === records[2] === records[3]
    @test sequence(records[1]) == "T"
    close(reader)

    # Test using new syntax
    filename = tempname()
    open(filename, "w") do io
        print(io, ">ABC\nUUGG\n>LLK\nAN\nPA\n\n")
    end
    Reader(NoopStream(open(filename))) do reader
        recs = collect(reader)
        @test map(identifier, recs) == ["ABC", "LLK"]
        @test map(sequence, recs) == ["UUGG", "ANPA"]
    end 
end

@testset "Reader edgecases" begin
    good = ">A\nA\n>B\nB"
    # do not accept > in sequence to detect missing newline
    # when two records are concatenated without newline
    bad = ">A\nA>B\nB"

    @test Reader(collect, IOBuffer(good)) isa Vector{Record}
    @test_throws Exception Reader(collect, IOBuffer(bad))
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
    writer = Writer(writer_buffer; width=0)
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

    # Test using new syntax
    filename = tempname()
    Writer(NoopStream(open(filename, "w"))) do writer
        write(writer, Record("some header", "UHAGC"))
        write(writer, Record("another_thing", "KJLM"))
    end
    @test read(filename, String) == ">some header\nUHAGC\n>another_thing\nKJLM\n"
end

# TranscodingStreams does not support flushing yet, so the FASTX implementation
# reaches into TS's internals. Hence we need to test it thoroughly
@testset "Writer flushing" begin
    strings = [
        ">hello there \nAACC\r\nTTGG",
        ">\n",
        ">someheader\n",
        "> tr wh [] ... ab **7\npwQ.---0l\r\n\r\naaaccc\r\n",
    ]
    target_strings = mktemp() do path, io
        map(strings) do string
            writer = Writer(open(path, "w"))
            write(writer, parse(Record, string))
            close(writer)
            open(io -> read(io, String), path)
        end
    end
    # Test that flushing at arbitrary points, then writing some more
    # works as expected
    for i in eachindex(strings)
        buf = IOBuffer()
        writer = Writer(buf)
        
        # First write some of the records and check that flush works
        for j in 1:i-1
            write(writer, parse(Record, strings[j]))
        end
        flush(writer)
        str = String(take!(copy(buf)))
        @test str == join(target_strings[1:i-1])

        # Then write the rest of them, and check the total results is as expected
        for j in i:lastindex(strings)
            write(writer, parse(Record, strings[j]))
        end
        flush(writer)
        str = String(take!(buf))
        @test str == join(target_strings)
    end
end

@testset "Writer width" begin
    header = "some data here"
    for width in (-10, 5, 25, 50)
        for seqlen in [width-1, width, 3*width, 3*width+3, 75, 200]
            seqlen < 1 && continue
            seq = join(rand('A':'Z', seqlen))
            record = Record(header, seq)
            buf = IOBuffer()
            writer = Writer(buf; width=width)
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
        writer = Writer(buf, width=0)
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

end # testset IO
