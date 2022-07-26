@testset "Valid specimens" begin
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
end

@testset "Invalid specimens" begin
    for specimen in list_invalid_specimens("FASTA")
        @test_throws Exception open(collect, Reader, path)
    end
end