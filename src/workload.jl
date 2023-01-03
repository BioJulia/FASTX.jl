using SnoopPrecompile: @precompile_setup, @precompile_all_calls

@precompile_setup begin
    fasta_path = joinpath(dirname(@__DIR__), "test", "test.fasta")
    fastq_path = joinpath(dirname(@__DIR__), "test", "test.fastq")
    fasta = read(fasta_path, String)
    fastq = read(fastq_path, String)

    @precompile_all_calls begin
        records = (
            parse(FASTA.Record, fasta),
            parse(FASTQ.Record, fastq)
        )
        for record in records
            identifier(record)
            description(record)
            seqsize(record)
            sequence(record)
            sequence(String, record)
            sequence(String, record)
            sequence(String, record)
        end

        # FASTQ specific
        record::FASTQ.Record = last(records)
        quality(record)
        collect(quality_scores(record))

        open(validate_fasta, fasta_path)
        open(validate_fastq, fastq_path)

        open(collect, FASTAReader, fasta_path)
        open(collect, FASTQReader, fastq_path)

        ind = open(faidx, fasta_path)
        rdr = FASTAReader(open(fasta_path); index=ind)
        extract(rdr, "abc", 2:3)
        seekrecord(rdr, "abc")
        close(rdr)
    end
end
