"""
    makelevels(levs::AbstractVector{<:AbstractString}[, inner=1, outer=1])

Create a `PooledVector` of the `levs` with repetitions `inner` and `outer`
"""
function makelevels(levs::AbstractVector{<:AbstractString}; inner=1, outer=1)
    repeat(PooledArray(levs), inner=inner, outer=outer)
end

function makelevels(nlev::Number, tag='S'; inner=1, outer=1)
    makelevels(string.(tag, lpad.(collect(1:nlev), ndigits(nlev), '0')), inner=inner,
        outer=outer)
end
