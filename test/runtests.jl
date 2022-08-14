using ReTest

#= If at some point we want to include in-line tests
using FASTX 
FASTX.retest() 
=#

include("TestFASTX.jl")
retest(TestFASTX)
include("TestFASTA.jl")
retest(TestFASTA)
include("TestFASTQ.jl")
retest(TestFASTQ)