using Documenter, FASTX

# Generate manual pages and examples.
include("generate.jl")

# Build documentation.

makedocs(
    format = :html,
    modules = [FASTX],
    sitename = "FASTX.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "FASTA formatted files" => "manual/generated/fasta.md",
            "FASTQ formatted files" => "manual/generated/fastq.md"
        ],
        "API Reference" => [
            "FASTA formatted files" => "reference/fasta.md",
            "FASTQ formatted files" => "reference/fastq.md"
        ]
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    julia = "1.0",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)