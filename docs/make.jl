using Documenter
using FASTX

# Build documentation.

makedocs(
    format = Documenter.HTML(),
    modules = [FASTX, FASTA, FASTQ],
    sitename = "FASTX.jl",
    doctest = false,
    strict = false,
    pages = [
        "FASTX" => "index.md",
        "FASTA" => "fasta.md",
        "FASTQ" => "fastq.md"
    ],
    authors = "Sabrina J. Ward, Jakob N. Nissen, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    deps = nothing,
    make = nothing
)