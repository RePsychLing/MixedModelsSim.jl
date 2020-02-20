using MixedModelsSim
using DataFrames, Tables
using Test

@testset "simdat_crossed both_btwn" begin
    # matching bwtn item/subj factors have different levels
    subj_btwn = Dict(:A => ["A", "B", "C"])
    item_btwn = Dict(:A => ["A", "B"])
    warn_str = "A has levels [\"A\", \"B\", \"C\"] for subj and [\"A\", \"B\"] for item"
    #@test_warn warn_str simdat_crossed(1, 1, subj_btwn = subj_btwn, item_btwn = item_btwn)


    for subj_n in [1, 5, 10, 20] 
        for item_n in [1, 5, 10, 20] 
            for b in [2, 3, 4, 5]
                btwn = nlevels(b, "A")
                subj_btwn = Dict(:A => btwn)
                item_btwn = Dict("A" => btwn)
                both_win = nothing

                dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn, item_btwn = item_btwn, both_win = both_win)

                @test nrow(dat) == subj_n*item_n*b
                @test ncol(dat) == 4
                @test names(dat) == [:subj, :A, :item, :dv]
                @test length(unique(dat.subj)) == subj_n*b
                @test length(unique(dat.item)) == item_n*b
            end
        end
    end

    for subj_n in [1, 5, 10, 20] 
        for item_n in [1, 5, 10, 20] 
            for b in [2, 3, 4, 5]
                btwn = nlevels(b, "A")
                # Dict keys can be strings or symbols, but have to be consistent
                subj_btwn = Dict("A" => btwn)
                item_btwn = Dict("A" => btwn)
                both_win = Dict(:time => ["pre", "post"])

                dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn, item_btwn = item_btwn, both_win = both_win)

                @test nrow(dat) == subj_n*item_n*b*2
                @test ncol(dat) == 5
                @test names(dat) == [:subj, :A, :item, :time, :dv]
                @test length(unique(dat.subj)) == subj_n*b
                @test length(unique(dat.item)) == item_n*b
            end
        end
    end

    # crossed between-item/subject factors
    subj_n = 1
    item_n = 1
    subj_btwn = Dict("age" => ["O", "Y"],
                     "A" => ["A1", "A2"],
                     "B" => ["B1", "B2"])
    item_btwn = Dict("A" => ["A1", "A2"],
                     "B" => ["B1", "B2"])
    both_win = nothing

    dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn, item_btwn = item_btwn, both_win = both_win)

    @test nrow(dat) == 8
    @test ncol(dat) == 6
    @test names(dat) == [:subj, :B, :A, :age, :item, :dv]
    @test dat.A == ["A1", "A1", "A2", "A2", "A1", "A1", "A2", "A2"]
    @test dat.B == ["B1", "B2", "B1", "B2", "B1", "B2", "B1", "B2"]
    @test dat.age == ["O", "O", "O", "O", "Y", "Y", "Y", "Y"]
    @test length(unique(dat.subj)) == 8
    @test length(unique(dat.item)) == 4

    # varying item_n and subj_n
    for subj_n in [1, 5, 10, 20]
        for item_n in [1, 5, 10, 20]
            subj_btwn = Dict(:age => ["O", "Y"],
                            :pet => ["cat", "dog"])
            item_btwn = Dict(:cond => ["A", "B"])
            both_win = Dict(:time => ["morning", "evening"])
            dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn, item_btwn = item_btwn, both_win = both_win)

            @test nrow(dat) == 16*subj_n*item_n
            @test length(unique(dat.subj)) == 4*subj_n
            @test length(unique(dat.item)) == 2*item_n
            @test names(dat) == [:subj, :age, :pet, :item, :cond, :time, :dv]
        end
    end

    # different numbers of levels for subj_btwn
    for n in 2:5
        subj_n = 10
        item_n = 10
        subj_btwn = Dict("n" => 1:n)
        dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn)

        @test nrow(dat) == n*subj_n*item_n
        @test length(unique(dat.subj)) == n*subj_n
        @test length(unique(dat.item)) == item_n
        @test names(dat) == [:subj, :n, :item, :dv]
    end

    # different numbers of levels for item_btwn
    for n in 2:5
        subj_n = 10
        item_n = 10
        item_btwn = Dict("n" => 1:n)
        dat = simdat_crossed(subj_n, item_n, item_btwn = item_btwn)

        @test nrow(dat) == n*subj_n*item_n
        @test length(unique(dat.subj)) == subj_n
        @test length(unique(dat.item)) == n*item_n
        @test names(dat) == [:subj, :item, :n, :dv]
    end

    # different numbers of levels for both_win
    for n in 2:5
        subj_n = 10
        item_n = 10
        both_win = Dict("n" => 1:n)
        dat = simdat_crossed(subj_n, item_n, both_win = both_win)

        @test nrow(dat) == n*subj_n*item_n
        @test length(unique(dat.subj)) == subj_n
        @test length(unique(dat.item)) == item_n
        @test names(dat) == [:subj, :item, :n, :dv]
    end

    for item_prefix in ["I", "item", "my item"]
        for subj_prefix in ["S", "subj", "my subj"]
            dat = simdat_crossed(2, 2, item_prefix = item_prefix, subj_prefix = subj_prefix)
            @test dat.subj == subj_prefix.*["1", "1", "2", "2"]
            @test dat.item == item_prefix.*["1", "2", "1", "2"]
        end
    end
end