using Documenter, FASTX

DocMeta.setdocmeta!(FASTX, :DocTestSetup, :(using FASTX, BioSequences); recursive=true)

makedocs(
    modules = [FASTX],
    format = Documenter.HTML(),
    sitename = "FASTX.jl",
    doctest = true,
    pages = [
        "Overview" => Any[
            "Overview" => "index.md",
            "Records" => "records.md",
            "File I/O" => "files.md",
        ],
        "FASTA" => "fasta.md",
        "FASTQ" => "fastq.md",
        "FAI" => "fai.md"
    ],
    authors = "Sabrina J. Ward, Jakob N. Nissen, The BioJulia Organisation and other contributors.",
    checkdocs = :all
)

deploydocs(
    repo = "github.com/BioJulia/FASTX.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)
