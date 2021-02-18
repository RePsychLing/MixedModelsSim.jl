"""
    simdat_crossed(RNG, subj_n, item_n; subj_btwn, item_btwn, both_win)
    simdat_crossed(subj_n, item_n; subj_btwn, item_btwn, both_win)

Return a `DataFrame` with a design specified by the:
* number of subjects (`subj_n`),
* number of items (`item_n`)
* between-subject factors (`subj_btwn`)
* between-item factors (`item_btwn`)
* within-subject/item factors (`both_win`)

If a factor is both between-subject and between-item,
put it in both `subj_btwn` and `item_btwn` with the same keys and the same levels.

Factors should be specified as Dicts in the following format:
```julia
Dict(
    :factor1_name => ["F1_level1", "F1_level2"],
    :factor2_name => ["F2_level1", "F2_level2", "F2_level3"]
)
```
Dict keys can be strings or symbols.
"""
function simdat_crossed(subj_n = 1, item_n = 1;
    subj_btwn = nothing, item_btwn = nothing, both_win = nothing,
    subj_prefix = "S", item_prefix = "I")


     simdat_crossed(Random.GLOBAL_RNG, subj_n, item_n;
                    subj_btwn = subj_btwn,
                    item_btwn = item_btwn,
                    both_win = both_win)
end

function simdat_crossed(rng::AbstractRNG, subj_n = 1, item_n = 1;
    subj_btwn = nothing, item_btwn = nothing, both_win = nothing,
    subj_prefix = "S", item_prefix = "I")

    # set up subj table
    if isnothing(subj_btwn)
        subj_vals = [nlevels(subj_n, subj_prefix)]
        sb_vars = []
    else
        sc = values(subj_btwn) |> collect
        sc = vcat([nlevels(subj_n)], sc)
        subj_prod = Iterators.product(sc...)
        subj_total_n = length(subj_prod)
        subj_vals = columntable(subj_prod) |> collect
        subj_vals[1] = nlevels(subj_total_n, subj_prefix) # rename subjects
        sb_vars = collect(keys(subj_btwn))
    end

    subj_names = Symbol.(vcat(["subj"], sb_vars))
    subj = NamedTuple{Tuple(subj_names)}(subj_vals)

    # set up item table
    if isnothing(item_btwn)
        item_vals = [nlevels(item_n, item_prefix)]
        ib_vars = []
    else
        ic = values(item_btwn) |> collect
        ic = vcat([nlevels(item_n)], ic)
        item_prod = Iterators.product(ic...)
        item_total_n = length(item_prod)
        item_vals = columntable(item_prod) |> collect
        item_vals[1] = nlevels(item_total_n, item_prefix) # rename items
        ib_vars = collect(keys(item_btwn))
    end

    item_names = Symbol.(vcat(["item"], ib_vars))
    item = NamedTuple{Tuple(item_names)}(item_vals)

    # set up within both table
    both_btwn_vars = Symbol.(intersect(subj_names, item_names))
    subj_df = DataFrame(subj)
    item_df = DataFrame(item)

    for bb in both_btwn_vars
        subj_levels = unique(subj[bb])
        item_levels = unique(item[bb])
        unique(item[bb])
        d = setdiff(subj_levels, item_levels)
        if length(d) > 0
            @warn(string(bb)*" has levels "*string(subj_levels)*" for subj and "*string(item_levels)*" for item")
        end
    end

    # cross the subject and item tables
    #design = factorproduct(subj, item) |> DataFrame
    if length(both_btwn_vars) == 0
        design = crossjoin(subj_df, item_df)
    else
        design = semijoin(subj_df, item_df; on = both_btwn_vars)
    end

    if (!isnothing(both_win))
        wc = values(both_win) |> collect
        win_prod = Iterators.product(wc...)
        win_vals = columntable(win_prod) |> collect
        win_names = collect(keys(both_win))
        win = NamedTuple{Tuple(Symbol.(win_names))}(win_vals)

        # cross the subject and item tables with any within factors
        #design = factorproduct(subj, item, win) |> DataFrame
        design = crossjoin(design, DataFrame(win))
    end

    # add random numbers as a DV
    design.dv = randn(rng, nrow(design))

    design

end

"""
    power_table(sim, alpha = 0.05)

Returns a `DataFrame` with two columns, `coefname` and `power`, with the proportion of
simulated p-values less than alpha, for `sim`, the output of `simulate_waldtests`.
"""
function power_table(sim, alpha = 0.05)
    pvals = DataFrame(columntable(sim).p)
    pvals = stack(pvals)
    pwr = combine(groupby(pvals, :variable), :value => x->mean(x.<alpha))
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
    rename!(widetab, ["iteration",  "coefname", "p",  "se",  "z",  "beta" ])
    sort!(widetab, [:iteration])
    select!(widetab, :iteration, :coefname, :beta, :se, :z, :p)
end
