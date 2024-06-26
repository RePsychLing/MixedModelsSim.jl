"""
    power_table(sim, alpha = 0.05)

Returns a row table with fields `coefname` and `power`  based on the proportion of
simulated p-values less than alpha, for `sim`, the output of `parametricbootstrap`.
"""
function power_table(sim::MixedModelBootstrap, alpha=0.05)
    # `se` represents the Bernoulli standard error of the number of trials below signifance using the (1 - alpha) for the probability.
    # q = 1 - alpha
    # se = sqrt(alpha * q / nsim)
    nsim = length(sim.objective)
    dd = Dict{Symbol,Int}()
    for row in sim.coefpvalues
        coef = row.coefname
        # make sure every even zero-power stuff is initialized
        val = get!(dd, coef, 0)
        if row.p < alpha
            dd[coef] = val + 1
        end
    end
    return [(; :coefname => string(key), :power => val / nsim) for (key, val) in dd]
end
