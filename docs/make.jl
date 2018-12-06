using Documenter, FASTX

# Generate manual pages and examples.
include("generate.jl")

# Build documentation.

makedocs(
    format = :html,
    modules = [FASTX, FASTX.FASTQ, FASTX.FASTA],
    sitename = "FASTX.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "FASTA formatted files" => "manual/generated/fasta.md",
            "FASTQ formatted files" => "manual/generated/fastq.md"
        ]
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    deps = nothing,
    make = nothing
)