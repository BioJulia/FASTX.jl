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
    @test quality(record) == "JJK"
    record = first(iterate(reader))
    @test identifier(record) == "B"
    @test description(record) == "B C"
    @test sequence(record) == "MNB"
    @test quality(record) == "LLL"
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
    records = [first(iterate(reader)) for i in 1:3]
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
    records = map(i -> parse(Record, i), [
        "@ABC DEF\nkjhmn\n+ABC DEF\njjjkk",
        "@pro [1-2](HLA=2);k=1\nttagga\n+\nabcdef",
    ])
    for record in records
        resize!(record.data, 512)
    end
    test_writer(records, r"@ABC DEF\nkjhmn\n\+ABC DEF\njjjkk\n@pro \[1-2\]\(HLA=2\);k=1\nttagga\n\+\nabcdef\n")

    # Exercise the IO with larger seqs to exceed the buffer
    target_buffer = IOBuffer()
    writer_buffer = IOBuffer()
    writer = Writer(writer_buffer)
    for i in 1:250
        name = join(rand('A':'z'), rand(20:30))
        name2 = join(rand('A':'z'), rand(30:40))
        descr = name * ' ' * name2
        seq = join(rand(('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'), rand(200:300)))
        qual_str = join(rand('A':'z', ncodeunits(seq)))
        write(target_buffer, '@', descr, '\n', seq, "\n+\n", qual_str, '\n')
        write(writer, Record(descr, seq, [Int8(i - OFFSET) for i in codeunits(qual_str)]))
    end
    flush(writer)
    writer_bytes = take!(writer_buffer)
    target_bytes = take!(target_buffer)
    close(writer)
    @test writer_bytes == target_bytes
end

# Rudamentary, see FASTA's tests
@testset "Writer flushing" begin
    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)
    target_strings = mktemp() do path, io
        map(records) do record
            writer = Writer(open(path, "w"))
            write(writer, record)
            close(writer)
            open(io -> read(io, String), path)
        end
    end

    for i in eachindex(records)
        buf = IOBuffer()
        writer = Writer(buf)
        
        # First write some of the records and check that flush works
        for j in 1:i-1
            write(writer, records[j])
        end
        flush(writer)
        str = String(take!(copy(buf)))
        @test str == join(target_strings[1:i-1])

        # Then write the rest of them, and check the total results is as expected
        for j in i:lastindex(records)
            write(writer,records[j])
        end
        flush(writer)
        str = String(take!(buf))
        @test str == join(target_strings)
    end
end

@testset "Writer optional second header" begin
    function iowrite(records, quality_header)
        buf = IOBuffer()
        writer = Writer(buf; quality_header=quality_header)
        for record in records
            write(writer, record)
        end
        flush(writer)
        str = String(take!(buf))
        close(writer)
        return str
    end

    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)

    @test iowrite(records, nothing) == join(map(string, records), '\n') * '\n'
    @test iowrite(records, true) == join(map(records) do record
        string(quality_header!(copy(record), true))
    end, '\n') * '\n'
    @test iowrite(records, false) == join(map(records) do record
        string(quality_header!(copy(record), false))
    end, '\n') * '\n'
end

@testset "Round trip" begin
    function iowrite(records)
        buf = IOBuffer()
        writer = Writer(buf)
        for record in records
            write(writer, record)
        end
        flush(writer)
        str = String(take!(buf))
        close(writer)
        return str
    end

    records = map(i -> parse(Record, i), TEST_RECORD_STRINGS)
    str = iowrite(records)
    records2 = Reader(collect, IOBuffer(str))
    str2 = iowrite(records2)
    @test str == str2
end