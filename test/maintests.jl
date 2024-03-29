# Common tests
@testset "FASTX" begin
    @testset "Copying to LongSequence" begin
        @test true
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
            # The next two have symbols that are in the FASTQ
            # accepted range A-z. So they can be parsed, but not
            # copied to sequence. If the FASTQ parsing changes,
            # just remove these from the test.
            "m^^jsommelig",        # Invalid
            "__m]]kvLMO",          # Invalid
        ]
        seqtypes = [
            LongDNA{4},
            LongDNA{2},
            LongRNA{4},
            LongRNA{2},
            LongAA
        ]
        for T in (FASTA.Record, FASTQ.Record)
            empty_record = if T == FASTA.Record
                record = T("name", "")
            else
                T("name", "", Int[])
            end 
            success = false
            seq = nothing
            for seqtype in seqtypes
                short_seq = seqtype()
                long_seq = seqtype(undef, 100)
                for str in strings
                    record = if T == FASTA.Record
                        record = T("name", str)
                    else
                        T("name", str, [75 for i in str])
                    end

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
    end

    @testset "Convert FASTQ to FASTA" begin
        for func in (FASTA.Record, i -> copy!(FASTA.Record(), i))
            rec = func(parse(FASTQ.Record, "@ta_g^ ha||;; \nTAGJKKm\n+\njjkkmmo"))
            @test description(rec) == "ta_g^ ha||;; "
            @test identifier(rec) == "ta_g^"
            @test sequence(rec) == "TAGJKKm"

            rec = func(parse(FASTQ.Record, "@\n\n+\n"))
            @test identifier(rec) == description(rec) == sequence(rec) == ""

            rec = func(parse(FASTQ.Record, "@mba M\npolA\n+mba M\nTAGA"))
            @test description(rec) == "mba M"
            @test identifier(rec) == "mba"
            @test sequence(rec) == "polA"
        end

        # Copyin conversion do no modify underlying record
        fq = parse(FASTQ.Record, "@ta_g^ ha||;; \nTAGJKKm\n+\njjkkmmo")
        fa = FASTA.Record(fq)
        fill!(fa.data, UInt8('a'))
        @test description(fq) == "ta_g^ ha||;; "
        @test identifier(fq) == "ta_g^"
        @test sequence(fq) == "TAGJKKm"
        @test quality(fq) == "jjkkmmo"
    end
end