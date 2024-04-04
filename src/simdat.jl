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

!!! note
    The number of subjects/items must divide the number of combinations of between subject/item factor levels.
    In other words, this function assumes a balanced design and will throw an error if that is not possible.
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
        sc = collect(values(subj_btwn))
        nlev = prod(length, sc)
        if mod(subj_n, nlev) != 0
            throw(ArgumentError("Number of subjects is not a multiple of the number of between-subject levels"))
        end
        subj_prod = Iterators.product(sc...)
        subj_vals = collect(columntable(subj_prod))
        subj_vals = vcat([nlevels(subj_n, subj_prefix)], repeat.(subj_vals, subj_n ÷ nlev))
        sb_vars = collect(keys(subj_btwn))
    end

    subj_names = vcat(["subj"], sb_vars)
    subj = (; (Symbol(k) => v for (k, v) in zip(subj_names, subj_vals))...)

    # set up item table
    if isnothing(item_btwn)
        item_vals = [nlevels(item_n, item_prefix)]
        ib_vars = []
    else
        ic = collect(values(item_btwn))
        nlev = prod(length, ic)
        if mod(item_n, nlev) != 0
            throw(ArgumentError("Number of items is not a multiple of the number of between-item levels"))
        end
        item_prod = Iterators.product(ic...)
        item_vals = collect(columntable(item_prod))
        item_vals = vcat([nlevels(item_n, item_prefix)], repeat.(item_vals, item_n ÷ nlev))
        ib_vars = collect(keys(item_btwn))
    end

    item_names = vcat(["item"], ib_vars)
    item = (; (Symbol(k) => v for (k, v) in zip(item_names, item_vals))...)

    # Check whether there are experimental factors which are both between-subject and between-item
    if isnothing(subj_btwn) || isnothing(item_btwn)
        both_between = []
    else
        both_between = intersect(keys(subj_btwn), keys(item_btwn))
    end

    # Case where there are not factors that are both within subject and within item
    if isnothing(both_win)
        # and there are no factors that are both between subject and between item
        if isempty(both_between)
            # cross the subject and item tables
            design = factorproduct(subj, item)
        else
            # make sure that each subject/item is only in one level of the between subject/between item factor
            design = [merge(x, y)
                      for x in rowtable(subj), y in rowtable(item)
                      if all(x[var] == y[var] for var in both_between)]
        end
    else
        # set up within both table
        wc = collect(values(both_win))
        win_prod = Iterators.product(wc...)
        win_vals = collect(columntable(win_prod))
        win_names = collect(keys(both_win))
        win = (; (Symbol(k) => v for (k, v) in zip(win_names, win_vals))...)

        if isempty(both_between)
            # cross the subject and item tables with any within factors
            design = factorproduct(subj, item, win)
        else
            # make sure that each subject/item is only in one level of the between subject/between item factor
            design = [merge(x, y, z)
                      for x in rowtable(subj), y in rowtable(item), z in rowtable(win)
                      if all(x[var] == y[var] for var in both_between)]
        end
    end

    dv = randn(rng, length(design))
    design = [merge(dd, (; :dv => vv)) for (dd, vv) in zip(design, dv)]

    return design
end
