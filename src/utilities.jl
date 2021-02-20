```@meta
DocTestSetup = quote
    using DataFrames, MixedModelsSim, Tables
end
```

"""
    factorproduct(facs...)

Return a `Vector{NamedTuple}` obtained by crossing `facs...`.

The arguments should be coerceable to a `Tables.RowTable` with `rowtable`.

The value is a `Tables.RowTable` and hence can be converted to a `DataFrame`.

# Example
```julia-repl
julia> DataFrame(factorproduct((item=nlevels(3,'I'),), (subj=nlevels(5), age=["Y","Y","Y","O","O"])))
15×3 DataFrame
│ Row │ item   │ subj   │ age    │
│     │ String │ String │ String │
├─────┼────────┼────────┼────────┤
│ 1   │ I1     │ S1     │ Y      │
│ 2   │ I2     │ S1     │ Y      │
│ 3   │ I3     │ S1     │ Y      │
│ 4   │ I1     │ S2     │ Y      │
│ 5   │ I2     │ S2     │ Y      │
│ 6   │ I3     │ S2     │ Y      │
│ 7   │ I1     │ S3     │ Y      │
│ 8   │ I2     │ S3     │ Y      │
│ 9   │ I3     │ S3     │ Y      │
│ 10  │ I1     │ S4     │ O      │
│ 11  │ I2     │ S4     │ O      │
│ 12  │ I3     │ S4     │ O      │
│ 13  │ I1     │ S5     │ O      │
│ 14  │ I2     │ S5     │ O      │
│ 15  │ I3     │ S5     │ O      │
```
"""
factorproduct(facs...) = factorproduct(Val(:rowtable), facs...)

function factorproduct(::Val{:rowtable}, facs...)
    vec([merge(s...) for s in Iterators.product(rowtable.(facs)...)])
end

# currently not more efficient, but reserving this use for later efficiency improvements
factorproduct(::Val{:columntable}, facs...) = columntable(factorproduct(facs...))

"""
    cyclicshift(v::AbstractVector, nrow)

Return an `eltype(v)` matrix of size `nrow` by `length(v)` with each column consisting
of `v` followed by a cyclic shift of `v` followed by ...

```@example
cyclicshift('a':'d', 8)
```
"""
function cyclicshift(v::AbstractVector, nrow)
    # The return value i used to counterbalance levels of conditions in [`withinitem`](@ref).
    vlen = length(v)
    [v[1 + (i + j) % vlen] for i in 0:(nrow - 1), j in 0:(vlen - 1)]
end

"""
    nlevels(nlev, tag='S')

Return a `Vector{String}` of `tag` followed by `1:nlev` left-padded with zeros

# Examples
```julia-repl
julia> show(nlevels(10))
["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10"]
```
"""
function nlevels(nlev, tag='S')
    string.(tag, lpad.(1:nlev, ndigits(nlev), '0'))
end

"""
    pooled!(df, cols::Type=Union{AbstractString,Missing})

Like `DataFrames.categorical!` but converting columns to `PooledArray`s

!!! warning
    This method is not type-specific in the first argument order to eliminate
    a dependency on `DataFrames.jl`. It nonetheless expects a `DataFrame` as
    its first argument
"""
function pooled!(df, cols::Type=Union{AbstractString,Missing})
    try
        for n in propertynames(df)
            v = df[!, n]
            if eltype(v) <: cols
                df[!, n] = PooledArray(v)
            end
        end
    catch e
        if e isa MethodError
            throw(ArgumentError("Something went wrong. Make sure you're using a DataFrame."))
        else
            rethrow(e)
        end
    end

    return df
end

"""
    flatlowertri(::LowerTriangular)

Returns the lower triangular flattened into 1D array in column-major order.
"""
function flatlowertri(l::LowerTriangular)
    rr, cc = size(l)
    return [l[i, j] for j = 1:cc for i = j:rr]
end


# """
#     withinitem(nitem, df)

# Return a `DataFrame` of `item` with `nitem` levels and balanced columns of conditions from `df`
# """
# function withinitem(nitem, df; tag = 'I')
#     nrow = size(df, 1)
#     q, r = divrem(nitem, nrow)
#     iszero(r) || throw(ArgumentError("nitem = $nitem is not a multiple of nrow = $nrow"))
#     value = df[vec(cyclicshift(1:nrow, nitem)), :]
#     value.item = repeat(PooledArray(nlevels(nitem, tag=tag)), outer=nrow)
#     value
# end
