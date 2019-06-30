using Documenter, FASTX

# Build documentation.

makedocs(
    format = Documenter.HTML(),
    modules = [FASTX, FASTX.FASTQ, FASTX.FASTA],
    sitename = "FASTX.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "FASTA formatted files" => "manual/fasta.md",
            "FASTQ formatted files" => "manual/fastq.md"
        ],
        "Library" => [
            "Public" => "lib/public.md"
        ]
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    deps = nothing,
    make = nothing
)