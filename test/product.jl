using MixedModelsSim
using Tables
using Test

@testset "factorproduct" begin
	items = (item = nlevels(50, 'I'), category = repeat(["ingroup", "outgroup"], inner=25))
	df = columntable(factorproduct(items, (subj = nlevels(100),)))
	@test length(df) == 3
	@test length(first(df)) == 5000
end
