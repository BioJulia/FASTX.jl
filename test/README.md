# Using ReTest for test development

1. If you don't already load it at startup, load [`Revise.jl`](https://github.com/timholy/Revise.jl)
2. Add the `test/` directory to your `LOAD_PATH`
3. Load ReTest and FASTX
4. Load the test module
5. Run tests early and often!

```julia-repl
❯ julia --project -q
julia> using Revise

julia> push!(LOAD_PATH, "./test")
4-element Vector{String}:
 "@"
 "@v#.#"
 "@stdlib"
 "./test"

julia> using ReTest, FASTX

julia> ReTest.load(FASTX)
Main.FASTXTests

julia> retest(FASTXTests)
                               Pass
Main.FASTXTests:
  FASTX                    |    975

Main.FASTXTests.TestFASTA:
  Record                   |     84
  IO                       |     52
  Index                    |    639
  Specimens                |     95
Main.FASTXTests.TestFASTA  |    870

Main.FASTXTests.TestFASTQ:
  Record                   |    161
                                    Pass    Fail   Total
  IO                            |     34       1      35
    Reader basics               |     14       1      15
    Writer basics               |      4               4
    Writer flushing             |     12              12
    Writer optional second header |      3               3
    Round trip                  |      1               1

Reader basics: Test Failed at /home/kevin/Repos/FASTX.jl/test/fastq/io.jl:4
  Expression: false

Main.FASTXTests.TestFASTQ       |    195       1     196
ERROR: Some tests did not pass: 34 passed, 1 failed, 0 errored, 0 broken.
```

Once you fix the test(s), re-running the command should only re-run the tests that have changed.

## Filtering

If you want to focus on a specific set of tests, you can pass one or more filters, eg

```julia-repl
julia> retest(FASTXTests, "basic")
                               Pass
Main.FASTXTests.TestFASTA:
  Record                   |     36
  IO                       |     18
Main.FASTXTests.TestFASTA  |     54

Main.FASTXTests.TestFASTQ:
  Record                   |     46
  IO                       |     18
Main.FASTXTests.TestFASTQ  |     64

Overall                    |    118

julia> retest(FASTXTests, "basic", "record")
                               Pass
Main.FASTXTests.TestFASTA:
  Record                   |     36

Main.FASTXTests.TestFASTQ:
  Record                   |     46

Overall                    |     82

julia> retest(FASTXTests, "specimens")
                               Pass
Main.FASTXTests.TestFASTA:
  Specimens                |     95

Main.FASTXTests.TestFASTQ:
  Specimens                |    115

Overall                    |    210
```