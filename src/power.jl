using MixedModels
using Random, StaticArrays, Tables

"""
    simulate_waldtests(rng::AbstractRNG, nsamp::Integer, m::LinearMixedModel;
        β = m.β, σ = m.σ, θ = m.θ, use_threads=false)
    simulate_waldtests(nsamp::Integer, m::LinearMixedModel;
        β = m.β, σ = m.σ, θ = m.θ, use_threads=false)

Return a Vector of βs, z- and p-values for each coefficient in a mixed model.
This is similar to the [`MixedModels.parametricbootstrap`](@ref), but returns
test statistics instead of estimates. This is useful for power analyses.

Perform `nsamp` parametric bootstrap replication fits of `m`,
returning a Vector of named tuples of the associated z- and p-values for the coefficents/

The default random number generator is `Random.GLOBAL_RNG`.

# Named Arguments

`β`, `σ`, and `θ` are the values of `m`'s parameters for simulating the responses.

`use_threads` determines whether threads are used to parallelize the computation.
Note that this functionality depends on the development version


```julia
using MixedModels, MixedModelsSim
using DataFrames, Gadfly, Random, StaticArrays, StatsBase, Tables

kb07 = MixedModels.dataset(:kb07);
form = @formula(rt_raw ~ 1 + spkr + prec + load + (1+spkr+prec+load|subj) + (1+spkr+prec+load|item));
cont = Dict(:spkr => HelmertCoding(),
            :prec => HelmertCoding(),
            :load => HelmertCoding())
fm1 = fit(MixedModel, form, kb07, contrasts=cont, REML=false);
zpmt = simulate_waldtests(MersenneTwister(42),10,fm1,use_threads=true);

mean(getindex.(columntable(zpmt).p, 1) .< 0.05)
mean(getindex.(columntable(zpmt).p, Symbol("(Intercept)")) .< 0.05)

mean(getindex.(columntable(zpmt).p, 2) .< 0.05)
mean(getindex.(columntable(zpmt).p, Symbol("spkr: old")) .< 0.05)

samples = DataFrame(columntable(zpmt).β)
shortestcovint(samples[!,Symbol("prec: maintain")])

plot(samples, layer(x=Symbol("load: yes")),
    Geom.density,
    Guide.xlabel("Bootstrap estimates of 'load: yes'"))

plot(stack(samples),
   layer(x=:value, color=:variable, xgroup=:variable),
   Geom.subplot_grid(Geom.density, free_x_axis=true),
   Guide.xlabel("Bootstrap estimates of β"))
```
"""
function simulate_waldtests(
    rng::AbstractRNG,
    n::Integer,
    morig::MixedModel{T};
    β = morig.β,
    σ = morig.σ,
    θ = morig.θ,
    use_threads = false,
) where {T}
    β = convert(Vector{T},β)
    σ = T(σ)
    θ = convert(Vector{T},θ)

    nβ, m = length(β), deepcopy(morig)
    # we need to do for in-place operations to work across threads
    m_threads = [m]
    coef_threads = [similar(β)]
    stderr_threads = [similar(β)]
    zstat_threads = [similar(β)]
    pval_threads = [similar(β)]

    if use_threads
        Threads.resize_nthreads!(m_threads)
        Threads.resize_nthreads!(coef_threads)
        Threads.resize_nthreads!(stderr_threads)
        Threads.resize_nthreads!(zstat_threads)
        Threads.resize_nthreads!(pval_threads)        
    end

    χsq = Chisq(1)

    cnms = Symbol.(coefnames(morig))

    tidynms = (:coefname,:beta,:se,:z,:p)
    results = NamedTuple{tidynms,Tuple{Symbol,T,T,T,T}}[]

    rnglock = ReentrantLock()
    pushlock = ReentrantLock()
    replicate(n, use_threads=use_threads) do
        mod = m_threads[Threads.threadid()]
        lock(rnglock)
        mod = simulate!(rng, mod, β = β, σ = σ, θ = θ)
        unlock(rnglock)
        refit!(mod)
        simβ = MixedModels.fixef!(coef_threads[Threads.threadid()], mod)
        simse = MixedModels.stderror!(stderr_threads[Threads.threadid()], mod)
        simz = map!(/, zstat_threads[Threads.threadid()], simβ, simse)
        simp = map!(z -> ccdf(χsq, abs2(z)), pval_threads[Threads.threadid()], simz)
        lock(pushlock)
        for t in zip(cnms, simβ, simse, simz, simp)
            push!(results, NamedTuple{tidynms}(t))
        end
        unlock(pushlock)
    end
    results
end

simulate_waldtests(
    n::Integer,
    morig::MixedModel{T};
    β = morig.β,
    σ = morig.σ,
    θ = morig.θ,
    use_threads = false,
) where {T} = simulate_waldtests(Random.GLOBAL_RNG, n, morig,
                                  β = β, σ = σ, θ = θ, use_threads = use_threads)
