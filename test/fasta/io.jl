
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
