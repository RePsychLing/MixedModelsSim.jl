# MixedModelsSim

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://RePsychLing.github.io/MixedModelsSim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://RePsychLing.github.io/MixedModelsSim.jl/dev)
![T1-url][T1-img]
![T2-url][T2-img]
[![Codecov](https://codecov.io/gh/RePsychLing/MixedModelsSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RePsychLing/MixedModelsSim.jl)

This package provides some utility functions for generating experimental designs, especially those with crossed factors.

## Installation

`MixedModelsSim` is not yet registered in the Julia package registry and must be installed by
```julia
julia> using Pkg
julia> Pkg.add(PackageSpec(url="https://github.com/RePsychLing/MixedModelsSim.jl"))
```
or, in the pkg REPL,
```julia
(@v1.4) pkg> add https://github.com/RePsychLing/MixedModelsSim.jl
```

## Purpose

This package provides functions to facilitate generating experimental designs, especially designs with crossed grouping factors such as "Subject" and "Item" in addition to experimental factors.  The experimental factors can be within-subject or within-item or between-subject and between-item.

This package uses structures from the [`Tables` package](https://github.com/JuliaData/Tables).  In particular, a data table can be viewed as a `rowtable`, which is a vector of `NamedTuple`s, or a `columntable` which is a `NamedTuple` of vectors (or something similar).

For those with experience in [`R`](https://www.r-project.org) just think of a `NamedTuple` as being like R's `list` type.  It's an ordered, named collection.

## Examples

To create a design with each of five subjects, three old and two young, tested on each of three items, first create the subject table
```julia
julia> using MixedModelsSim, DataFrames, Tables
julia> subject = (subj = ["S1","S2","S3","S4","S5"], age=["O","O","O","Y","Y"]);
julia> typeof(subject)
NamedTuple{(:subj, :age),Tuple{Array{String,1},Array{String,1}}}

julia> rowtable(subject)
5-element Array{NamedTuple{(:subj, :age),Tuple{String,String}},1}:
 (subj = "S1", age = "O")
 (subj = "S2", age = "O")
 (subj = "S3", age = "O")
 (subj = "S4", age = "Y")
 (subj = "S5", age = "Y")
 ```
 then create the design as the product of an item table (defined inline here) and the `subject` table
 ```julia
julia> design = factorproduct((item = ["I1","I2","I3"],), subject)
15-element Array{NamedTuple{(:item, :subj, :age),Tuple{String,String,String}},1}:
 (item = "I1", subj = "S1", age = "O")
 (item = "I2", subj = "S1", age = "O")
 (item = "I3", subj = "S1", age = "O")
 (item = "I1", subj = "S2", age = "O")
 (item = "I2", subj = "S2", age = "O")
 (item = "I3", subj = "S2", age = "O")
 (item = "I1", subj = "S3", age = "O")
 (item = "I2", subj = "S3", age = "O")
 (item = "I3", subj = "S3", age = "O")
 (item = "I1", subj = "S4", age = "Y")
 (item = "I2", subj = "S4", age = "Y")
 (item = "I3", subj = "S4", age = "Y")
 (item = "I1", subj = "S5", age = "Y")
 (item = "I2", subj = "S5", age = "Y")
 (item = "I3", subj = "S5", age = "Y")
```

The design can be converted to a `DataFrame` and the strings pooled to save storage.
```julia
julia> design |> DataFrame |> pooled!
15×3 DataFrame
│ Row │ item   │ subj   │ age    │
│     │ String │ String │ String │
├─────┼────────┼────────┼────────┤
│ 1   │ I1     │ S1     │ O      │
│ 2   │ I2     │ S1     │ O      │
│ 3   │ I3     │ S1     │ O      │
│ 4   │ I1     │ S2     │ O      │
│ 5   │ I2     │ S2     │ O      │
│ 6   │ I3     │ S2     │ O      │
│ 7   │ I1     │ S3     │ O      │
│ 8   │ I2     │ S3     │ O      │
│ 9   │ I3     │ S3     │ O      │
│ 10  │ I1     │ S4     │ Y      │
│ 11  │ I2     │ S4     │ Y      │
│ 12  │ I3     │ S4     │ Y      │
│ 13  │ I1     │ S5     │ Y      │
│ 14  │ I2     │ S5     │ Y      │
│ 15  │ I3     │ S5     │ Y      │

julia> describe(ans)
3×8 DataFrame
│ Row │ variable │ mean    │ min    │ median  │ max    │ nunique │ nmissing │ eltype   │
│     │ Symbol   │ Nothing │ String │ Nothing │ String │ Int64   │ Nothing  │ DataType │
├─────┼──────────┼─────────┼────────┼─────────┼────────┼─────────┼──────────┼──────────┤
│ 1   │ item     │         │ I1     │         │ I3     │ 3       │          │ String   │
│ 2   │ subj     │         │ S1     │         │ S5     │ 5       │          │ String   │
│ 3   │ age      │         │ O      │         │ Y      │ 2       │          │ String   │
```

## Background on tables and tuples

In Julia tuples are created by listing the contents, surrounded by parentheses and separated by commas.

```julia
julia> Tables.istable(subject)
true

julia> Tables.schema(subject)
Tables.Schema:
 :subj  String
 :age   String

 julia> DataFrame(subject)
 5×2 DataFrame
 │ Row │ subj   │ age    │
 │     │ String │ String │
 ├─────┼────────┼────────┤
 │ 1   │ S1     │ Y      │
 │ 2   │ S2     │ Y      │
 │ 3   │ S3     │ Y      │
 │ 4   │ S4     │ O      │
 │ 5   │ S5     │ O      │
 ```
### The curiously trailing comma

To distinguish creating a named tuple of length 1 from an assignment with parentheses around it a comma is required after the first named element.  To create an item table with only three item identifiers, the expression must be written
```julia
julia> items = (item = ["I1", "I2", "I3"],)
(item = ["I1", "I2", "I3"],)
```
with that curiously trailing comma.  In general, trailing commas are allowed in the creation of tuples or in argument lists but in this case the trailing comma is mandatory.

## Generating factors with n levels

The `nlevels` utility function can be used to generate a vector of length `n` with a given tag.  For example, the vector of subject levels can be generated as
```julia
julia> show(nlevels(5, 'S'))
["S1", "S2", "S3", "S4", "S5"]
```
The default tag is `S` so this sequence could be generated more simply as
```julia
julia> show(nlevels(5))
["S1", "S2", "S3", "S4", "S5"]
```

[T1-img]: https://github.com/JuliaStats/MixedModels.jl/workflows/Tier1/badge.svg
[T1-url]: https://github.com/JuliaStats/MixedModels.jl/actions?workflow=Tier1

[T2-img]: https://github.com/JuliaStats/MixedModels.jl/workflows/Tier2/badge.svg
[T2-url]: https://github.com/JuliaStats/MixedModels.jl/actions?workflow=Tier2
