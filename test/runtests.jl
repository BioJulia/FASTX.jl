using ReTest

#= If at some point we want to include in-line tests
using FASTX 
FASTX.retest() 
=#

include("FASTXTests.jl")
retest(FASTXTests)
