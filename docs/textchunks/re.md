## Assemble the Random Effects

The hard part in simulating right now is specifying the random effects.
We're working on making this bit easier, but you need to specify the variance-covariance matrix of the random effects. You can see what this
looks like:

```@example Main
vc = VarCorr(m0)
```

```@example Main
vc.σρ
```

In matrix form:

```@example Main
m0.λ
```

```@example Main
m0.λ[1]
```

```@example Main
m0.λ[2]
```

Note that the entries in `m0.λ` are stored in the same order as the random effects are listed in `VarCorr(m0)`.  The off-diagonal elements are covariances and correspond to the correlations (the ρ's).
The on-diagonal elements are just the standard deviations (the σ's).
For this example, we'll just assume all the correlations and thus the covariances are 0 in order to make things simple.
Then we only have to worry about the diagonal elements.

Let's assume that the variability
- between items
    - in the intercept is 1.3 times the residual variability
    - in age is 0.35 times the residual variability
    - in context is 0.75 times the residual variability
- between subjects
    - in the intercept is 1.5 times the residual variability
    - in frequency is 0.5 times the residual variability
    - in context is 0.75 times the residual variability

Then we'll have the following λ matrices:

```@example Main
λitem = LowerTriangular(diagm([1.3, 0.35, 0.75]))
```

```@example Main
λsubj = LowerTriangular(diagm([1.5, 0.5, 0.75]))
```

For the actual optimization process, MixedModels.jl actually uses a flattened version of these stored in the θ vector. (More precisely, these λ matrices are derived from the θ vector.)

```@example Main
isapprox(m0.θ,  [flatlowertri(m0.λ[1]); flatlowertri(m0.λ[2])])
```

With that in mind, we can assemble our θ vector for simulation:

```@example Main
# make sure these are in the same order as in the model summary!
θ = [flatlowertri(λitem); flatlowertri(λsubj)]
```

We can check that we got these right by installing these parameter values into the model:
