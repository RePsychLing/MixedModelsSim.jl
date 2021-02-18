"""
    replicate(f::Function, n::Integer; use_threads=false)

Return a vector of the values of `n` calls to `f()` - used in simulations where the value of `f` is stochastic.

Note that if `f()` is not thread-safe or depends on a non thread-safe RNG,
    then you must set `use_threads=false`. Also note that ordering of replications
    is not guaranteed when `use_threads=true`, although the replications are not
    otherwise affected for thread-safe `f()`.
"""
function replicate(f::Function, n::Integer; use_threads=false)
    if use_threads
        # no macro version yet: https://github.com/timholy/ProgressMeter.jl/issues/143
        p = Progress(n)
        # get the type
        rr = f()
        next!(p)
        # pre-allocate
        results = [rr for _ in Base.OneTo(n)]
        Threads.@threads for idx = 2:n
            results[idx] = f()
            next!(p)
        end
    else
        results = @showprogress [f() for _ in Base.OneTo(n)]
    end
    results
end