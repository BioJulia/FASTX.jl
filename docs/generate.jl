import Literate

# Generate the manual pages

MANUALDIR = joinpath(@__DIR__, "src", "manual")
MANUALOUT = joinpath(@__DIR__, "src", "manual", "generated")

for manpage in readdir(MANUALDIR)
    endswith(manpage, ".jl") || continue
    input = abspath(joinpath(MANUALDIR, manpage))
    Literate.markdown(input, MANUALOUT)
end