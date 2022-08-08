using Documenter, FASTX

makedocs(
    format = Documenter.HTML(),
    sitename = "FASTX.jl",
    doctest = true,
    pages = [
        "FASTX" => "index.md",
        "FASTA" => "fasta.md",
        "FASTQ" => "fastq.md",
        "FAI" => "fai.md"
    ],
    authors = "Sabrina J. Ward, Jakob N. Nissen, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)
