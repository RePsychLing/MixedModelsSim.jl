# Understanding tasks

Feel free to ask questions. 

Do not assume that requested tasks are possible: feel free to inform me that a requested task is not possible in the given form or that there are only tedious workarounds.

# Julia development

For Julia code, use the YAS codestyle (https://github.com/jrevels/YASGuide). 
Explicit `return` statements are required and the use of `import` is forbidden.
Place all `using` statements in the top-level module file or in `set_up_tests.jl`.

Always invoke Julia with `--startup-file=no` unless explicitly instructed otherwise.

You may invoke a specific Julia version with `+VERSION`, e.g. `+1.10` or `+1.12`. This argument must come immediately after `julia` and before any other flags. 

Don't forget to specify the local project with the `--project=.`.

## Writing tests

When checking test coverage, you can use LocalCoverage.jl, which writes coverage to `coverage/lcov.info`.

When running individual tests in Julia, you need to load the test environment. This can be done with `using TestEnv; TestEnv.activate()`. The entire testsuite can be invoked with `using Pkg; Pkg.test()` without activating the test environment beforehand, but you will still need to activate the local project (e.g. with the `--project=.` argument to Julia).

Do not use `@test_warn`. When testing warnings, use `@test_logs` with the appropriate logging level, e.g. `:warn`. For tests where a warning is issued, use the `@suppressor` macro from Suppressor.jl to hide the warning during testing.

All shared test state — `using` statements, shared constants, and shared helper functions — belongs in `set_up_tests.jl`. Individual test files should not define their own constants or helpers if they are used across files.

When creating temporary directories in tests, use the `do`-block form: `mktempdir() do tempdir ... end`. This ensures cleanup even if the test errors.
