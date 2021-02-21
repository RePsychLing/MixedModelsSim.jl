using MixedModelsSim
using DataFrames, Tables
using Test

@testset "simdat_crossed" begin
    subj_btwn = Dict(:age => ["O", "Y"],
                     :pet => ["cat", "dog"])
    item_btwn = Dict(:cond => ["A", "B"])
    both_win = Dict(:time => ["morning", "evening"])

    for n in [4, 8]
        subj_n = n
        item_n = n
        # since we don't inspect the dv here, we don't worry about specifying the RNG
        # this has the advantage of testing the default/implicit-RNG method
        dat = simdat_crossed(subj_n, item_n,
                             subj_btwn = subj_btwn,
                             item_btwn = item_btwn,
                             both_win = both_win,
                             subj_prefix="SS",
                             item_prefix="LL")

        @test Tables.isrowtable(dat)

        @test first(dat).subj[1:2] == "SS"
        @test first(dat).item[1:2] == "LL"

        @test length(dat) == 2 * subj_n * item_n
        @test length(unique(getproperty.(dat, :subj))) == subj_n
        @test length(unique(getproperty.(dat, :item))) == item_n
        @test Set(Tables.columnnames(first(dat))) == Set([:subj, :pet, :age, :item, :cond, :time, :dv])
    end

    @test_throws ArgumentError simdat_crossed(3, 1, subj_btwn = subj_btwn)
    @test_throws ArgumentError simdat_crossed(1, 3, item_btwn = item_btwn)

    # the nothingness
    dat = simdat_crossed(1, 1)
    @test Tables.isrowtable(dat)
    @test Set(Tables.columnnames(first(dat))) == Set([:subj, :item, :dv])
end
