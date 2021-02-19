using MixedModelsSim
using DataFrames, Tables
using Test

@testset "simdat_crossed" begin
    for n in [1, 5, 10, 20]
        subj_n = n
        item_n = n
        subj_btwn = Dict(:age => ["O", "Y"],
                         :pet => ["cat", "dog"])
        item_btwn = Dict(:cond => ["A", "B"])
        both_win = Dict(:time => ["morning", "evening"])
        # since we don't inspect the dv here, we don't worry about specifying the RNG
        # this has the advantage of testing the default/implicit-RNG method
        dat = simdat_crossed(subj_n, item_n,
                             subj_btwn = subj_btwn,
                             item_btwn = item_btwn,
                             both_win = both_win)

        @test Tables.isrowtable(dat)

        @test length(dat) == 16 * subj_n * item_n
        @test length(unique(getproperty.(dat, :subj))) == 4 * subj_n
        @test length(unique(getproperty.(dat, :item))) == 2 * item_n
        @test Set(Tables.columnnames(first(dat))) == Set([:subj, :pet, :age, :item, :cond, :time, :dv])
    end

    # different numbers of levels for subj_btwn
    for n in 2:5
        subj_n = 10
        item_n = 10
        subj_btwn = Dict(:n => 1:n)
        dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn)

        @test length(dat) == n * subj_n * item_n
        @test length(unique(getproperty.(dat, :subj))) == n * subj_n
        @test length(unique(getproperty.(dat, :item))) == item_n
        @test Set(Tables.columnnames(first(dat))) == Set([:subj, :n, :item, :dv])
    end

    # different numbers of levels for item_btwn
    for n in 2:5
        subj_n = 10
        item_n = 10
        item_btwn = Dict(:n => 1:n)
        dat = simdat_crossed(subj_n, item_n, item_btwn = item_btwn)

        @test length(dat) == n * subj_n * item_n
        @test length(unique(getproperty.(dat, :subj))) == subj_n
        @test length(unique(getproperty.(dat, :item))) == n * item_n
        @test Set(Tables.columnnames(first(dat))) == Set([:subj, :n, :item, :dv])
    end

    # different numbers of levels for both_win
    for n in 2:5
        subj_n = 10
        item_n = 10
        both_win = Dict(:n => 1:n)
        dat = simdat_crossed(subj_n, item_n, both_win = both_win)

        @test length(dat) == n * subj_n * item_n
        @test length(unique(getproperty.(dat, :subj))) == subj_n
        @test length(unique(getproperty.(dat, :item))) == item_n
        @test Set(Tables.columnnames(first(dat))) == Set([:subj, :n, :item, :dv])
    end
end
