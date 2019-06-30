# Only run coverage from linux nightly build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "1.1" || exit()

using Coverage

cd(joinpath(@__DIR__, "..")) do
    Codecov.submit(Codecov.process_folder())
end