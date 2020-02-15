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
factorproduct(facs...) = vec([merge(s...) for s in Iterators.product(rowtable.(facs)...)])

"""
    cyclicshift(v::AbstractVector, nrow)

Return an `eltype(v)` matrix of size `nrow` by `length(v)` with each column consisting
of `v` followed by a cyclic shift of `v` followed by ...

The return value is used to counterbalance levels of conditions in [`withinitem`](@ref).
```@example
cyclicshift('a':'d', 8)
```
"""
function cyclicshift(v::AbstractVector, nrow)
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
nlevels(nlev, tag='S') = string.(tag, lpad.(1:nlev, ndigits(nlev), '0'))

"""
    pooled!(df::DataFrame, cols::Type=Union{AbstractString,Missing})

Like `DataFrames.categorical!` but converting columns to `PooledArray`s
"""
function pooled!(df::DataFrame, cols::Type=Union{AbstractString,Missing})
    for (n,v) in eachcol(df, true)
        if eltype(v) <: cols
            df[!, n] = PooledArray(v)
        end
    end
    df
end

"""
    withinitem(nitem, df)

Return a `DataFrame` of `item` with `nitem` levels and balanced columns of conditions from `df`
"""
function withinitem(nitem, df; tag = 'I')
    nrow = size(df, 1)
    q, r = divrem(nitem, nrow)
    iszero(r) || throw(ArgumentError("nitem = $nitem is not a multiple of nrow = $nrow"))
    value = df[vec(cyclicshift(1:nrow, nitem)), :]
    value.item = repeat(PooledArray(nlevels(nitem, tag=tag)), outer=nrow)
    value
end

"""
    power_table(sim, 0.05)

Returns a `DataFrame` with two columns, `coefname` and `power`, with the proportion of 
simulated p-values less than alpha, for `sim`, the output of `simulate_waldtests`.
"""

function power_table(sim, alpha = 0.05)
    pvals = DataFrame(columntable(sim).p)
    pvals = stack(pvals) 
    pwr = by(pvals, :variable, :value => x->mean(x.<alpha) )
    rename!(pwr, ["coefname", "power"])
end

"""
    sim_to_df(sim)

Returns a `DataFrame` with 6 columns: `iteration`, `coefname`, `beta`, `se`, `z`, `p`. 
Rows are all the coefficients for each iteration of `sim`, the output of `simulate_waldtests`. 
`iteration` is not guaranteed to be the same across runs of `simulate_waldtests` with the same seed, 
even though the samples will be.
"""

function sim_to_df(sims)
    tab = DataFrame()
    for (i, sim) in enumerate(sims)
        df = DataFrame(sim)
        df[!, :iteration] .= i
        df[!, :var] .= string.(collect(keys(sim)))
        append!(tab, df)
    end
    longtab = stack(tab, 1:(ncol(tab)-2), variable_name = :coefname)
    widetab = unstack(longtab, :var, :value)
    rename!(widetab, ["coefname", "iteration",  "p",  "se",  "z",  "beta" ])
    df_ordered = widetab[[:iteration, :coefname, :beta, :se, :z, :p]]
    sort!(df_ordered, [:iteration])
end
