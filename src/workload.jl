using SnoopPrecompile

@precompile_all_calls begin
    fasta = ">abc def\nTAG\nTA"
    fastq = "@ABC DEF\nTAGC\n+ABC DEF\nJJJJ"
    records = (
        parse(FASTA.Record, fasta),
        parse(FASTQ.Record, fastq)
    )
    for record in records
        identifier(record)
        description(record)
        seqsize(record)
        sequence(record)
        sequence(BioSequences.LongDNA{2}, record)
        sequence(BioSequences.LongDNA{4}, record)
        sequence(BioSequences.LongAA, record)
    end

    # FASTQ specific
    record::FASTQ.Record = last(records)
    quality(record)
    collect(quality_scores(record))

    validate_fasta(IOBuffer(fasta))
    validate_fastq(IOBuffer(fastq))

    collect(FASTAReader(IOBuffer(fasta)))
    collect(FASTQReader(IOBuffer(fastq)))

    ind = faidx(IOBuffer(fasta))
    rdr = FASTAReader(IOBuffer(fasta); index=ind)
    extract(rdr, "abc", 2:3)
    seekrecord(rdr, "abc")
end
