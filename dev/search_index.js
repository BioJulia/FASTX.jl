var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#FASTX-1",
    "page": "Home",
    "title": "FASTX",
    "category": "section",
    "text": "(Image: Latest Release) (Image: MIT license)  (Image: Stable documentation) (Image: Latest documentation) (Image: Lifecycle) (Image: Chat on Discord)"
},

{
    "location": "#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "FASTX provides a Reader for I/O for the FASTA, text based sequence format. "
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "FASTX is bundled into the Bio.jl package, so you may not need to install this package explicitly. However, if you do, you can install BioSequences from the Julia REPL:using Pkg\nadd(\"FASTX\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "#Testing-1",
    "page": "Home",
    "title": "Testing",
    "category": "section",
    "text": "FASTA is tested against Julia 1.X on Linux, OS X, and Windows.Latest release Latest build status\n(Image: ) (Image: ) (Image: )(Image: )"
},

{
    "location": "#Contributing-1",
    "page": "Home",
    "title": "Contributing",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.Take a look at the CONTRIBUTING file provided with every BioJulia package package for detailed contributor and maintainer guidelines."
},

{
    "location": "#Financial-contributions-1",
    "page": "Home",
    "title": "Financial contributions",
    "category": "section",
    "text": "We also welcome financial contributions in full transparency on our open collective. Anyone can file an expense. If the expense makes sense for the development of the community, it will be \"merged\" in the ledger of our open collective by the core contributors and the person who filed the expense will be reimbursed."
},

{
    "location": "#Backers-and-Sponsors-1",
    "page": "Home",
    "title": "Backers & Sponsors",
    "category": "section",
    "text": "Thank you to all our backers and sponsors!Love our work and community? Become a backer.(Image: backers)Does your company use BioJulia? Help keep BioJulia feature rich and healthy by sponsoring the project Your logo will show up here with a link to your website.(Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: )"
},

{
    "location": "#Questions?-1",
    "page": "Home",
    "title": "Questions?",
    "category": "section",
    "text": "If you have a question about contributing or using BioJulia software, come on over and chat to us on Discord, or you can try the Bio category of the Julia discourse site."
},

{
    "location": "manual/generated/fasta/#",
    "page": "FASTA formatted files",
    "title": "FASTA formatted files",
    "category": "page",
    "text": "EditURL = \"https://github.com/BioJulia/FASTX.jl/blob/master/docs/src/manual/fasta.jl\"CurrentModule = FASTX"
},

{
    "location": "manual/generated/fasta/#IO-FASTA-formatted-files-1",
    "page": "FASTA formatted files",
    "title": "IO - FASTA formatted files",
    "category": "section",
    "text": "FASTA is a text-based file format for representing biological sequences. A FASTA file stores a list of sequence records with name, description, and sequence.The template of a sequence record is:>{name} {description}?\n{sequence}Here is an example of a chromosomal sequence:>chrI chromosome 1\nCCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC\nCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.Reader",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.Reader",
    "category": "type",
    "text": "FASTA.Reader(input::IO; index = nothing)\n\nCreate a data reader of the FASTA file format.\n\nArguments\n\ninput: data source\nindex=nothing: filepath to a random access index (currently fai is supported)\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.Writer",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.Writer",
    "category": "type",
    "text": "FASTA.Writer(output::IO; width=70)\n\nCreate a data writer of the FASTA file format.\n\nArguments\n\noutput: data sink\nwidth=70: wrapping width of sequence characters\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.Record",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.Record",
    "category": "type",
    "text": "FASTA.Record()\n\nCreate an unfilled FASTA record.\n\n\n\n\n\nFASTA.Record(data::Vector{UInt8})\n\nCreate a FASTA record object from data.\n\nThis function verifies and indexes fields for accessors. Note that the ownership of data is transferred to a new record object.\n\n\n\n\n\nFASTA.Record(str::AbstractString)\n\nCreate a FASTA record object from str.\n\nThis function verifies and indexes fields for accessors.\n\n\n\n\n\nFASTA.Record(identifier, sequence)\n\nCreate a FASTA record object from identifier and sequence.\n\n\n\n\n\nFASTA.Record(identifier, description, sequence)\n\nCreate a FASTA record object from identifier, description and sequence.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.hasidentifier",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.hasidentifier",
    "category": "function",
    "text": "hasidentifier(record::Record)\n\nChecks whether or not the record has an identifier.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.identifier",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.identifier",
    "category": "function",
    "text": "identifier(record::Record)::String\n\nGet the sequence identifier of record.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.hasdescription",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.hasdescription",
    "category": "function",
    "text": "hasdescription(record::Record)\n\nChecks whether or not the record has a description.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.description",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.description",
    "category": "function",
    "text": "description(record::Record)::String\n\nGet the description of record.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.hassequence",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.hassequence",
    "category": "function",
    "text": "hassequence(record::Record)\n\nChecks whether or not a sequence record contains a sequence.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#FASTX.FASTA.sequence-Tuple{FASTX.FASTA.Record,UnitRange{Int64}}",
    "page": "FASTA formatted files",
    "title": "FASTX.FASTA.sequence",
    "category": "method",
    "text": "sequence(record::Record, [part::UnitRange{Int}])\n\nGet the sequence of record.\n\nThis function infers the sequence type from the data. When it is wrong or unreliable, use sequence(::Type{S}, record::Record).  If part argument is given, it returns the specified part of the sequence.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fasta/#Readers-and-Writers-1",
    "page": "FASTA formatted files",
    "title": "Readers and Writers",
    "category": "section",
    "text": "The reader and writer for FASTA formatted files, are found within the BioSequences.FASTA submodule.FASTA.Reader\nFASTA.WriterThey can be created with IOStreams:using FASTX\n\nr = FASTA.Reader(open(\"my-reads.fasta\", \"r\"))\nw = FASTA.Writer(open(\"my-out.fasta\", \"w\"))Alternatively, Base.open is overloaded with a method for conveinience:r = open(FASTA.Reader, \"my-reads.fasta\")\nw = open(FASTA.Writer, \"my-out.fasta\")Usually sequence records will be read sequentially from a file by iteration.reader = open(FASTA.Reader, \"my-reads.fasta\")\nfor record in reader\n    # Do something\nend\nclose(reader)But if the FASTA file has an auxiliary index file formatted in fai, the reader supports random access to FASTA records, which would be useful when accessing specific parts of a huge genome sequence:reader = open(FASTA.Reader, \"sacCer.fa\", index = \"sacCer.fa.fai\")\nchrIV = reader[\"chrIV\"]  # directly read sequences called chrIV.\nclose(reader)Reading in a sequence from a FASTA formatted file will give you a variable of type FASTA.Record.FASTA.RecordVarious getters and setters are available for FASTA.Records:FASTA.hasidentifier\nFASTA.identifier\nFASTA.hasdescription\nFASTA.description\nFASTA.hassequence\nFASTA.sequence(record::FASTA.Record, [part::UnitRange{Int}])To write a BioSequence to FASTA file, you first have to create a FASTA.Record:using BioSequences\nx = dna\"aaaaatttttcccccggggg\"\nrec = FASTA.Record(\"MySeq\", x)\nw = open(FASTA.Writer, \"my-out.fasta\")\nwrite(w, rec)\nclose(w)As always with julia IO types, remember to close your file readers and writer after you are finished.Using open with a do-block can help ensure you close a stream after you are finished.open(FASTA.Reader, \"my-reads.fasta\") do reader\n    for record in reader\n        # Do something\n    end\nendThis page was generated using Literate.jl."
},

{
    "location": "manual/generated/fastq/#",
    "page": "FASTQ formatted files",
    "title": "FASTQ formatted files",
    "category": "page",
    "text": "EditURL = \"https://github.com/BioJulia/FASTX.jl/blob/master/docs/src/manual/fastq.jl\"CurrentModule = FASTX"
},

{
    "location": "manual/generated/fastq/#IO-FASTQ-formatted-files-1",
    "page": "FASTQ formatted files",
    "title": "IO - FASTQ formatted files",
    "category": "section",
    "text": "FASTQ is a text-based file format for representing DNA sequences along with qualities for each base. A FASTQ file stores a list of sequence records in the following format:@{name} {description}?\n{sequence}\n+\n{qualities}Here is an example of one record from a FASTQ file:@FSRRS4401BE7HA\ntcagTTAAGATGGGAT\n+\n###EEEEEEEEE##E#"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.Reader",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.Reader",
    "category": "type",
    "text": "FASTQ.Reader(input::IO; fill_ambiguous=nothing)\n\nCreate a data reader of the FASTQ file format.\n\nArguments\n\ninput: data source\nfill_ambiguous=nothing: fill ambiguous symbols with the given symbol\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.Writer",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.Writer",
    "category": "type",
    "text": "FASTQ.Writer(output::IO; quality_header=false)\n\nCreate a data writer of the FASTQ file format.\n\nArguments\n\noutput: data sink\nquality_header=false: output the title line at the third line just after \'+\'\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.Record",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.Record",
    "category": "type",
    "text": "FASTQ.Record()\n\nCreate an unfilled FASTQ record.\n\n\n\n\n\nFASTQ.Record(data::Vector{UInt8})\n\nCreate a FASTQ record object from data.\n\nThis function verifies and indexes fields for accessors. Note that the ownership of data is transferred to a new record object.\n\n\n\n\n\nFASTQ.Record(str::AbstractString)\n\nCreate a FASTQ record object from str.\n\nThis function verifies and indexes fields for accessors.\n\n\n\n\n\nFASTQ.Record(identifier, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, sequence and quality.\n\n\n\n\n\nFASTQ.Record(identifier, description, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, description, sequence and quality.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.hasidentifier",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.hasidentifier",
    "category": "function",
    "text": "hasidentifier(record::Record)\n\nChecks whether or not the record has an identifier.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.identifier",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.identifier",
    "category": "function",
    "text": "identifier(record::Record)::String\n\nGet the sequence identifier of record.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.hasdescription",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.hasdescription",
    "category": "function",
    "text": "hasdescription(record::Record)\n\nChecks whether or not the record has a description.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.description",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.description",
    "category": "function",
    "text": "description(record::Record)::String\n\nGet the description of record.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.hassequence",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.hassequence",
    "category": "function",
    "text": "hassequence(record::Record)\n\nChecks whether or not a sequence record contains a sequence.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.sequence-Tuple{FASTX.FASTQ.Record,UnitRange{Int64}}",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.sequence",
    "category": "method",
    "text": "sequence(record::Record, [part::UnitRange{Int}])::BioSequences.DNASequence\n\nGet the sequence of record.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.hasquality",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.hasquality",
    "category": "function",
    "text": "hasquality(record::Record)\n\nCheck whether the given FASTQ record has a quality string.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.quality",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.quality",
    "category": "function",
    "text": "quality(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}\n\nGet the base quality of record.\n\n\n\n\n\nquality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}\n\nGet the base quality of record by decoding with encoding_name.\n\nThe encoding_name can be either :sanger, :solexa, :illumina13, :illumina15, or :illumina18.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#FASTX.FASTQ.Record-Tuple{AbstractString,Union{Nothing, AbstractString},Any,Array{T,1} where T}",
    "page": "FASTQ formatted files",
    "title": "FASTX.FASTQ.Record",
    "category": "method",
    "text": "FASTQ.Record(identifier, description, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, description, sequence and quality.\n\n\n\n\n\n"
},

{
    "location": "manual/generated/fastq/#Readers-and-Writers-1",
    "page": "FASTQ formatted files",
    "title": "Readers and Writers",
    "category": "section",
    "text": "The reader and writer for FASTQ formatted files, are found within the BioSequences.FASTQ module.FASTQ.Reader\nFASTQ.WriterThey can be created with IOStreams:using FASTX\n\nr = FASTQ.Reader(open(\"../my-reads.fastq\", \"r\"))\nw = FASTQ.Writer(open(\"my-output.fastq\", \"w\"))Alternatively, Base.open is overloaded with a method for conveinience:r = open(FASTQ.Reader, \"../my-reads.fastq\")\nw = open(FASTQ.Writer, \"my-out.fastq\")Note that FASTQ.Reader does not support line-wraps within sequence and quality. Usually sequence records will be read sequentially from a file by iteration.reader = open(FASQ.Reader, \"../my-reads.fastq\")\nfor record in reader\n    # Do something\nend\nclose(reader)Reading in a record from a FASTQ formatted file will give you a variable of type FASTQ.Record.FASTQ.RecordVarious getters and setters are available for FASTQ.Records:FASTQ.hasidentifier\nFASTQ.identifier\nFASTQ.hasdescription\nFASTQ.description\nFASTQ.hassequence\nFASTQ.sequence(record::FASTQ.Record, [part::UnitRange{Int}])\nFASTQ.hasquality\nFASTQ.qualityTo write a BioSequence to FASTQ file, you first have to create a FASTQ.Record:FASTQ.Record(identifier::AbstractString, description::Union{AbstractString,Nothing}, sequence, quality::Vector; offset=33)As always with julia IO types, remember to close your file readers and writer after you are finished.Using open with a do-block can help ensure you close a stream after you are finished.open(FASTQ.Reader, \"../my-reads.fastq\") do reader\n    for record in reader\n        # Do something\n    end\nend#-This page was generated using Literate.jl."
},

{
    "location": "reference/fasta/#",
    "page": "FASTA formatted files",
    "title": "FASTA formatted files",
    "category": "page",
    "text": "CurrentModule = FASTX"
},

{
    "location": "reference/fasta/#FASTA-1",
    "page": "FASTA formatted files",
    "title": "FASTA",
    "category": "section",
    "text": ""
},

{
    "location": "reference/fasta/#FASTA-Reader-1",
    "page": "FASTA formatted files",
    "title": "FASTA Reader",
    "category": "section",
    "text": "FASTA.Reader"
},

{
    "location": "reference/fasta/#FASTA-Record-1",
    "page": "FASTA formatted files",
    "title": "FASTA Record",
    "category": "section",
    "text": "FASTA.Record\nFASTA.identifier\nFASTA.hasidentifier\nFASTA.description\nFASTA.hasdescription\nFASTA.sequence"
},

{
    "location": "reference/fasta/#FASTA-Writer-1",
    "page": "FASTA formatted files",
    "title": "FASTA Writer",
    "category": "section",
    "text": "FASTA.Writer"
},

{
    "location": "reference/fastq/#",
    "page": "FASTQ formatted files",
    "title": "FASTQ formatted files",
    "category": "page",
    "text": "CurrentModule = FASTX"
},

{
    "location": "reference/fastq/#FASTQ-1",
    "page": "FASTQ formatted files",
    "title": "FASTQ",
    "category": "section",
    "text": ""
},

{
    "location": "reference/fastq/#FASTQ-Reader-1",
    "page": "FASTQ formatted files",
    "title": "FASTQ Reader",
    "category": "section",
    "text": "FASTQ.Reader"
},

{
    "location": "reference/fastq/#FASTQ-Record-1",
    "page": "FASTQ formatted files",
    "title": "FASTQ Record",
    "category": "section",
    "text": "FASTQ.Record\nFASTQ.identifier\nFASTQ.hasidentifier\nFASTQ.description\nFASTQ.hasdescription\nFASTQ.sequence\nFASTQ.hassequence\nFASTQ.quality\nFASTQ.hasquality\nFASTQ.QualityEncoding\nFASTQ.SANGER_QUAL_ENCODING\nFASTQ.SOLEXA_QUAL_ENCODING\nFASTQ.ILLUMINA13_QUAL_ENCODING\nFASTQ.ILLUMINA15_QUAL_ENCODING\nFASTQ.ILLUMINA18_QUAL_ENCODING "
},

{
    "location": "reference/fastq/#FASTQ-Writer-1",
    "page": "FASTQ formatted files",
    "title": "FASTQ Writer",
    "category": "section",
    "text": "FASTQ.Writer"
},

]}
