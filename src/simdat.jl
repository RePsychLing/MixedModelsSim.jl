"""
    simdat_crossed([RNG], subj_n, item_n;
                   subj_btwn=nothing, item_btwn=nothing, both_win=nothing,
                   subj_prefix="S", item_prefix="I")

Return a row table with a design specified by the:
* number of subjects (`subj_n`),
* number of items (`item_n`)
* between-subject factors (`subj_btwn`)
* between-item factors (`item_btwn`)
* within-subject/item factors (`both_win`)

If a factor is both between-subject and between-item, put it in both
`subj_btwn` and `item_btwn` with the same keys and the same levels.

Factors should be specified as dictionaries in the following format:
```julia
Dict(
    :factor1_name => ["F1_level1", "F1_level2"],
    :factor2_name => ["F2_level1", "F2_level2", "F2_level3"]
)
```

In addition to design, the rowtable contains a field `dv` pre-populated
with N(0,1) noise as a basis for further simulating a design.
"""
simdat_crossed(args...; kwargs...) = simdat_crossed(Random.GLOBAL_RNG, args...; kwargs...)

function simdat_crossed(rng::AbstractRNG, subj_n=1, item_n=1;
    subj_btwn=nothing, item_btwn=nothing, both_win=nothing,
    subj_prefix="S", item_prefix="I")

    # set up subj table
    if isnothing(subj_btwn)
        subj_vals = [nlevels(subj_n, subj_prefix)]
        sb_vars = []
    else
        sc = values(subj_btwn) |> collect
        nlev = prod(length, sc)
        if mod(subj_n, nlev) != 0
            throw(ArgumentError("Number of subjects is not a multiple of the number of between-subject levels"))
        end
        subj_prod = Iterators.product(sc...)
        subj_vals = columntable(subj_prod) |> collect
        subj_vals = vcat([nlevels(subj_n, subj_prefix)], repeat.(subj_vals, subj_n ÷ nlev))
        sb_vars = collect(keys(subj_btwn))
    end

    subj_names = vcat(["subj"], sb_vars)
    subj = (; (Symbol(k) => v for (k,v) in zip(subj_names, subj_vals))...)

    # set up item table
    if isnothing(item_btwn)
        item_vals = [nlevels(item_n, item_prefix)]
        ib_vars = []
    else
        ic = values(item_btwn) |> collect
        nlev = prod(length, ic)
        if mod(item_n, nlev) != 0
            throw(ArgumentError("Number of items is not a multiple of the number of between-subject levels"))
        end
        item_prod = Iterators.product(ic...)
        item_vals = columntable(item_prod) |> collect
        item_vals = vcat([nlevels(item_n, item_prefix)], repeat.(item_vals, item_n ÷ nlev))
        ib_vars = collect(keys(item_btwn))
    end

    item_names = vcat(["item"], ib_vars)
    item = (; (Symbol(k) => v for (k,v) in zip(item_names, item_vals))...)

    # set up within both table
    if (isnothing(both_win))
        # cross the subject and item tables
        design = factorproduct(subj, item)
    else
        wc = values(both_win) |> collect
        win_prod = Iterators.product(wc...)
        win_vals = columntable(win_prod) |> collect
        win_names = collect(keys(both_win))
        win = (; (Symbol(k) => v for (k,v) in zip(win_names, win_vals))...)

        # cross the subject and item tables with any within factors
        design = factorproduct(subj, item, win)
    end

    dv = randn(rng, length(design))
    design = [ merge(dd, (; :dv => vv) ) for (dd, vv) in zip(design, dv) ]

    design
end
