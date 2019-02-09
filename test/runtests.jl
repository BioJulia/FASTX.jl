using Test
using FASTX
using FormatSpecimens
import BioGenerics
import BioGenerics.Testing: intempdir
import BioCore
import BioCore.Testing.bio_fmt_specimens
import BioSequences:
    @dna_str,
    @rna_str,
    @aa_str,
    DNASequence,
    BioSequence,
    DNAAlphabet,
    DNA_N,
    DNA_A,
    DNA_G

@testset "FASTA" begin
    @testset "Record" begin
        record = FASTA.Record()
        @test !BioGenerics.isfilled(record)
        @test_throws ArgumentError FASTA.identifier(record)

        record = FASTA.Record(">foo\nACGT\n")
        @test BioGenerics.isfilled(record)
        @test BioGenerics.hasseqname(record)
        @test FASTA.hasidentifier(record)
        @test BioGenerics.seqname(record) == FASTA.identifier(record) == "foo"
        @test !FASTA.hasdescription(record)
        @test_throws BioGenerics.Exceptions.MissingFieldException FASTA.description(record)
        @test BioGenerics.hassequence(record)
        @test FASTA.hassequence(record)
        @test FASTA.sequence(record) == dna"ACGT"
        @test FASTA.sequence(record, 2:3) == dna"CG"
        @test FASTA.sequence(String, record) == "ACGT"
        @test FASTA.sequence(String, record, 2:3) == "CG"
        @test record == FASTA.Record(">foo\nACGT\n")

        record = FASTA.Record("""
        >CYS1_DICDI fragment
        SCWSFSTTGNVEGQHFISQNKL
        VSLSEQNLVDCDHECMEYEGE
        """)
        @test BioGenerics.isfilled(record)
        @test FASTA.identifier(record) == "CYS1_DICDI"
        @test FASTA.description(record) == "fragment"
        @test FASTA.sequence(record) == aa"SCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGE"
        @test FASTA.sequence(record, 10:15) == aa"NVEGQH"
    end

    output = IOBuffer()
    writer = FASTA.Writer(output, 5)
    write(writer, FASTA.Record("seq1", dna"TTA"))
    write(writer, FASTA.Record("seq2", "some description", dna"ACGTNN"))
    flush(writer)
    @test String(take!(output)) == """
    >seq1
    TTA
    >seq2 some description
    ACGTN
    N
    """

    reader = FASTA.Reader(IOBuffer(
    """
    >seqA some description
    QIKDLLVSSSTDLDTTLKMK
    ILELPFASGDLSM
    >seqB
    VLMALGMTDLFIPSANLTG*
    """))
    record = FASTA.Record()
    @test read!(reader, record) === record
    @test FASTA.identifier(record) == "seqA"
    @test FASTA.description(record) == "some description"
    @test FASTA.sequence(record) == aa"QIKDLLVSSSTDLDTTLKMKILELPFASGDLSM"
    @test read!(reader, record) === record
    @test FASTA.identifier(record) == "seqB"
    @test !FASTA.hasdescription(record)
    @test FASTA.sequence(record) == aa"VLMALGMTDLFIPSANLTG*"

    function test_fasta_parse(filename, valid)
        # Reading from a stream
        stream = open(FASTA.Reader, filename)
        @test eltype(stream) == FASTA.Record
        if valid
            for seqrec in stream end
            @test true  # no error
            @test close(stream) === nothing
        else
            @test_throws Exception begin
                for seqrec in stream end
            end
            return
        end

        # in-place parsing
        stream = open(FASTA.Reader, filename)
        entry = eltype(stream)()
        while !eof(stream)
            read!(stream, entry)
        end

        # Check round trip
        output = IOBuffer()
        writer = FASTA.Writer(output, width = 60)
        outputB = IOBuffer()
        writerB = FASTA.Writer(outputB, width = -1)
        expected_entries = Any[]
        for seqrec in open(FASTA.Reader, filename)
            write(writer, seqrec)
            write(writerB, seqrec)
            push!(expected_entries, seqrec)
        end
        flush(writer)
        flush(writerB)

        seekstart(output)
        seekstart(outputB)
        read_entries = FASTA.Record[]
        read_entriesB = FASTA.Record[]
        for seqrec in FASTA.Reader(output)
            push!(read_entries, seqrec)
        end
        for seqrec in FASTA.Reader(outputB)
            push!(read_entriesB, seqrec)
        end
        @test expected_entries == read_entries
        @test expected_entries == read_entriesB
    end

    function valid_specimen_filter(specimen)
        tags = specimen["tags"]
        valid = get(specimen, "valid", true)
        return valid && !occursin("comments", tags)
    end
    function invalid_specimen_filter(specimen)
        tags = specimen["tags"]
        valid = get(specimen, "valid", true)
        return !valid && !occursin("comments", tags)
    end
    
    fasta_folder = path_of_format("FASTA")     
    valid_specimens = list_valid_specimens("FASTA") do specimen
        !hastag(specimen, "comments")
    end
    invalid_specimens = list_invalid_specimens("FASTA")
    for specimen in valid_specimens
        test_fasta_parse(joinpath(fasta_folder, filename(specimen)), true)
    end
    for specimen in invalid_specimens
        test_fasta_parse(joinpath(fasta_folder, filename(specimen)), false)
    end

    @testset "Faidx" begin
        fastastr = """
        >chr1
        CCACACCACACCCACACACC
        >chr2
        ATGCATGCATGCAT
        GCATGCATGCATGC
        >chr3
        AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
        TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
        >chr4
        TACTT
        """
        # generated with `samtools faidx`
        faistr = """
        chr1	20	6	20	21
        chr2	28	33	14	15
        chr3	100	69	50	51
        chr4	5	177	5	6
        """
        mktempdir() do dir
            filepath = joinpath(dir, "test.fa")
            write(filepath, fastastr)
            write(filepath * ".fai", faistr)
            open(FASTA.Reader, filepath, index=filepath * ".fai") do reader
                chr3 = reader["chr3"]
                @test FASTA.identifier(chr3) == "chr3"
                @test FASTA.sequence(chr3) == dna"""
                AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
                TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
                """

                chr2 = reader["chr2"]
                @test FASTA.identifier(chr2) == "chr2"
                @test FASTA.sequence(chr2) == dna"""
                ATGCATGCATGCAT
                GCATGCATGCATGC
                """

                chr4 = reader["chr4"]
                @test FASTA.identifier(chr4) == "chr4"
                @test FASTA.sequence(chr4) == dna"""
                TACTT
                """

                chr1 = reader["chr1"]
                @test FASTA.identifier(chr1) == "chr1"
                @test FASTA.sequence(chr1) == dna"""
                CCACACCACACCCACACACC
                """

                @test_throws ArgumentError reader["chr5"]
            end
        end

        # invalid index
        @test_throws ArgumentError FASTA.Reader(IOBuffer(fastastr), index=π)
    end

    @testset "append" begin
        intempdir() do
            filepath = "test.fa"
            writer = open(FASTA.Writer, filepath)
            write(writer, FASTA.Record("seq1", dna"AAA"))
            close(writer)
            writer = open(FASTA.Writer, filepath, append=true)
            write(writer, FASTA.Record("seq2", dna"CCC"))
            close(writer)
            seqs = open(collect, FASTA.Reader, filepath)
            @test length(seqs) == 2
        end
    end
end

@testset "FASTQ" begin
    @test isa(FASTQ.Record("1", dna"AA", UInt8[10, 11]), FASTQ.Record)
    @test isa(FASTQ.Record("1", "desc.", dna"AA", UInt8[10, 11]), FASTQ.Record)
    @test_throws ArgumentError FASTQ.Record("1", dna"AA", UInt8[10])

    output = IOBuffer()
    writer = FASTQ.Writer(output)
    write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
    write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
    flush(writer)
    @test String(take!(output)) == """
    @1
    AN
    +
    ,:
    @2 high quality
    TGA
    +
    IJN
    """

    output = IOBuffer()
    writer = FASTQ.Writer(output, quality_header=true)
    write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
    write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
    flush(writer)
    @test String(take!(output)) == """
    @1
    AN
    +1
    ,:
    @2 high quality
    TGA
    +2 high quality
    IJN
    """

    @testset "Record" begin
        record = FASTQ.Record()
        @test !BioCore.isfilled(record)

        record = FASTQ.Record("""
        @SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
        AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
        +SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
        @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
        """)
        @test BioCore.isfilled(record)
        @test FASTQ.hasidentifier(record) == BioCore.hasseqname(record) == true
        @test FASTQ.identifier(record) == BioCore.seqname(record) == "SRR1238088.1.1"
        @test FASTQ.hasdescription(record)
        @test FASTQ.description(record) == "HWI-ST499:111:D0G94ACXX:1:1101:1173:2105"
        @test FASTQ.hassequence(record) == BioCore.hassequence(record) == true
        @test FASTQ.sequence(DNASequence, record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.sequence(record) == BioCore.sequence(record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.sequence(String, record) == "AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.hasquality(record)
        @test FASTQ.quality(record) == b"@BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ" .- 33

        record = FASTQ.Record("""
        @SRR1238088.1.1
        AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
        +
        @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
        """)
        @test BioCore.isfilled(record)
        @test !FASTQ.hasdescription(record)
    end

    function test_records(rs1, rs2)
        if length(rs1) != length(rs2)
            return false
        end
        for (r1, r2) in zip(rs1, rs2)
            if FASTQ.identifier(r1) != FASTQ.identifier(r2) ||
               FASTQ.sequence(r1)   != FASTQ.sequence(r2)   ||
               FASTQ.quality(r1)    != FASTQ.quality(r2)
                return false
            end
        end
        return true
    end

    function test_fastq_parse(filename, valid)
        # Reading from a reader
        reader = open(FASTQ.Reader, filename)
        @test eltype(reader) == FASTQ.Record
        if valid
            for record in reader end
            @test true  # no error
            @test close(reader) === nothing
        else
            @test_throws Exception begin
                for record in reader end
            end
            return
        end

        # in-place parsing
        reader = open(FASTQ.Reader, filename)
        record = eltype(reader)()
        try
            while true
                read!(reader, record)
            end
        catch ex
            close(reader)
            if !isa(ex, EOFError)
                rethrow()
            end
        end

        # Check round trip
        output = IOBuffer()
        writer = FASTQ.Writer(output)
        expected_entries = FASTQ.Record[]
        for record in open(FASTQ.Reader, filename)
            write(writer, record)
            push!(expected_entries, record)
        end
        flush(writer)

        seekstart(output)
        read_entries = FASTQ.Record[]
        for record in FASTQ.Reader(output)
            push!(read_entries, record)
        end

        return test_records(expected_entries, read_entries)
    end

    function valid_specimen_filter(specimen)
        tags = split(get(specimen, "tags", ""))
        if any(t ∈ tags for t in ["gaps", "rna", "comments", "linewrap"])
            return false
        end
        valid = get(specimen, "valid", true)
        return valid
    end
    function invalid_specimen_filter(specimen)
        tags = split(get(specimen, "tags", ""))
        if any(t ∈ tags for t in ["gaps", "rna", "comments", "linewrap"])
            return false
        end
        valid = get(specimen, "valid", true)
        return !valid
    end
         
    valid_specimens = bio_fmt_specimens("FASTQ", valid_specimen_filter)
    invalid_specimens = bio_fmt_specimens("FASTQ", invalid_specimen_filter)
    #println("testing valid fastqs")
    for specimen in valid_specimens
        #println(specimen)
        test_fastq_parse(specimen, true)
    end
    #println("testing invalid fastqs")
    for specimen in invalid_specimens
        #println(specimen)
        test_fastq_parse(specimen, false)
    end

    @testset "invalid quality encoding" begin
        # Sanger full range (note escape characters before '$' and '\')
        record = FASTQ.Record("""
        @FAKE0001 Original version has PHRED scores from 0 to 93 inclusive (in that order)
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
        +
        !"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        """)

        # the range is not enough in these encodings
        for encoding in (:solexa, :illumina13, :illumina15)
            @test_throws ErrorException FASTQ.quality(record, encoding)
        end

        # the range is enough in these encodings
        for encoding in (:sanger, :illumina18)
            @test FASTQ.quality(record, encoding) == collect(0:93)
        end
    end

    @testset "fill ambiguous nucleotides" begin
        input = IOBuffer("""
        @seq1
        ACGTNRacgtnr
        +
        BBBB##AAAA##
        """)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=nothing))) == dna"ACGTNRACGTNR"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_A)))   == dna"ACGTAAACGTAA"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_G)))   == dna"ACGTGGACGTGG"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_N)))   == dna"ACGTNNACGTNN"
        seekstart(input)
        @test FASTQ.sequence(BioSequence{DNAAlphabet{2}}, first(FASTQ.Reader(input, fill_ambiguous=DNA_A))) == dna"ACGTAAACGTAA"
    end
end

@testset "Quality scores" begin
    @testset "Decoding base quality scores" begin
        function test_decode(encoding, values, expected)
            result = Array{Int8}(undef, length(expected))
            FASTQ.decode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected
        end

        test_decode(FASTQ.SANGER_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])

        test_decode(FASTQ.SOLEXA_QUAL_ENCODING,
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[-5, 2, 3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA13_QUAL_ENCODING,
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA15_QUAL_ENCODING,
                    UInt8['C', 'D', 'E', 'h', '~'],
                    Int8[3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA18_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])
    end

    @testset "Encoding base quality scores" begin
        function test_encode(encoding, values, expected)
            result = Array{UInt8}(undef, length(expected))
            FASTQ.encode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected
        end

        test_encode(FASTQ.SANGER_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])

        test_encode(FASTQ.SOLEXA_QUAL_ENCODING,
                    Int8[-5, 2, 3, 4, 5, 40, 62],
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA13_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 62],
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA15_QUAL_ENCODING,
                    Int8[3, 4, 5, 40, 62],
                    UInt8['C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA18_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])
    end
end

