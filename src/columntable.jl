"""
    nlevels(nlev, tag='S')

Return a `Vector{String}` of `tag` followed by `1:nlev` padded with zeros
"""
nlevels(nlev, tag='S') = string.(tag, lpad.(1:nlev, ndigits(nlev), '0'))

"""
    crossedDataFrame(namedtup)

Return a `DataFrame` obtained by crossing factors given in `namedtup`.

`namedtup` should be a named tuple of vectors of strings giving the levels.
"""
function crossedDataFrame(namedtup)
    nms = keys(namedtup)
    rowtbl = NamedTuple{nms, NTuple{length(namedtup), String}}[]
    for pr in Iterators.product(values(namedtup)...)
        push!(rowtbl, NamedTuple{nms}(pr))
    end
    ctbl = Tables.columntable(rowtbl)
    DataFrame(NamedTuple{keys(ctbl)}(PooledArray.(values(ctbl))))
end

function itemsubjdf(nitem, nsubj, expfactors=())
    df = crossedDataFrame((item=nlevels(nitem, 'I'), subj=nlevels(nsubj)))
    if !isempty(expfactors)
        exptbl = crossedDataFrame(expfactors)
        ncond = size(exptbl, 1)
        if !iszero(rem(nitem, ncond))
            @warn "number of items, $nitem, is not a multiple of number of conditions, $ncond"
        end
        if !iszero(rem(nsubj, ncond))
            @warn "number of subjects, $nsubj, is not a multiple of number of conditions, $ncond"
        end
        nrow = size(df, 1)
        q, r = divrem(nrow, ncond)
        if !iszero(r)
            throw(ArgumentError("number of rows, $nrow, is not a multiple of number of conditions, $ncond"))
        end
        df = hcat(df, repeat(exptbl, q))
    end
    df
end

function cyclicshift!(v::AbstractVector)
    tmp = first(v)
    @inbounds for i in 2:length(v)
        v[i-1] = v[i]
    end
    v[end] = tmp
    v
end
