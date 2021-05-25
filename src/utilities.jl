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
    nlevstbl(nm::Symbol, n, vars::Pair{Symbol, Vector{String}}...)

Return a `Tables.columntable` with a `nm` column as a `PooledArray` with `n` levels.

If any `vars` pairs are given they are expanded to columns representing characteristics
of the `nm` column.  In experimental design terminology, if say `nm` is `:item` then these
represent between-item experimental factors.

The `nm` column is generated as `nlevels(n, uppercase(first(string(nm))))`

# Examples
```julia-repl
julia> nlevstbl(:item, 10)
(item = ["I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I10"],)

julia> nlevstbl(:item, 9, :level => ["low", "medium", "high"])
(item = ["I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9"], level = ["low", "medium", "high", "low", "medium", "high", "low", "medium", "high"])
```
"""
function nlevstbl(nm::Symbol, n::Integer, vars::Pair{Symbol, Vector{String}}...)
    nms = [nm]
    vals = [PooledArray(nlevels(n, uppercase(first(string(nm)))), signed=true, compress=true)]
    inner = 1
    for var in vars
        levs = last(var)
        nlev = length(levs)
        rept = inner * nlev
        q, r = divrem(n, rept)
        iszero(r) || throw(ArgumentError("n = $n is not a multiple of repetition block size $rept"))
        push!(vals, PooledArray(repeat(levs, inner=inner, outer=q), signed=true, compress=true))
        push!(nms, first(var))
        inner = rept
    end
    NamedTuple{(nms...,)}((vals...,))
end

"""
    pooled!(df, cols::Type=Union{AbstractString,Missing})

Like `DataFrames.categorical!` but converting columns to `PooledArray`s

!!! warning
    This method is not type-specific in the first argument, in order to eliminate
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


"""
    update!(m::MixedModel; θ)

Update the mixed model to use θ as its new parameter vector.

!!! note
    This is a convenience function for installing a particular parameter vector
    and the resulting model fit. It does not actually perform any type of
    optimization.

!!! note
    For GLMMs, this only sets θ and not β, even for `fast=false` fits.
"""
update!(m::LinearMixedModel; θ) = updateL!(MixedModels.setθ!(m, θ))
update!(m::GeneralizedLinearMixedModel; θ) = pirls!(MixedModels.setθ!(m, θ), false)
# arguably type piracy, but we're all the same developers....



"""
    update!(m::MixedModel, re...)

Update the mixed model to use the random-effects covariance matrices.

The `re` can be created using [`create_re`](@ref).

They should be specified in the order specified in `VarCorr(m)`.

Details
========
The `re` used as the λ fields of the model's `ReTerm`s and should be specified
as the lower Cholesky factor of covariance matrices.
"""
function update!(m::MixedModel, re...)
    θ = vcat((flatlowertri(rr) for rr in re)...)
    update!(m; θ=θ)
end

"""
    create_re(sigmas...; corrmat=Matrix{Float64}(I, length(sigmas), length(sigmas))

Create the covariance factor for a random effect from the standard deviations and correlation matrix.

The `sigmas` should be specified in the same order as the random slopes in the
output of `VarCorr(m)`.

The correlation matrix defaults to the identiy matrix, i.e. no correlation between
random effects.

!!! note
    The return value is the lower Cholesky factor of the covariance matrix, which is what
    [`update!`](@ref) requires.
"""
function create_re(sigmas...; corrmat=nothing)
    ss = Diagonal([sigmas...])

    isnothing(corrmat) ? LowerTriangular(ss) : ss * cholesky(Symmetric(corrmat, :L)).L
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
