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

@testset "simdat_crossed between-subjects between-items - simple case" begin
    # stimulate data for a case in which a factor is both between-subject and between-item
    conditions = Dict(:cond => ["A", "B"])
    cond_n = length(conditions[:cond])
    subj_n = 2
    item_n = 4

    data = DataFrame(simdat_crossed(subj_n, item_n,
                                    subj_btwn = conditions,
                                    item_btwn = conditions))

    @test nrow(data) == subj_n * (item_n/cond_n)

    # Test whether each subject is only in one of the conditions
    for s in nlevels(subj_n,"S")
        @test length(unique(data[isequal.(data.subj, s), :cond])) == 1
    end

    # Test whether each item is only in one of the conditions
    for i in nlevels(item_n,"I")
        @test length(unique(data[isequal.(data.item, i), :cond])) == 1
    end

end;

@testset "simdat_crossed between-subjects between-items - complex case" begin
    # stimulate data for a case in which a factor is both between-subject and between-item
    both_btwn = Dict(:cond1 => ["A", "B"],
                     :cond2 => ["C", "D", "E"])
    subj_btwn = merge(Dict(:age => ["O", "Y"]), both_btwn)                    
    item_btwn = merge(Dict(:pet => ["cat", "dog"]), both_btwn)
    both_win = Dict(:time => ["morning", "evening"])

    subj_n = 12
    item_n = 12

    data = DataFrame(simdat_crossed(subj_n, item_n,
        subj_btwn = subj_btwn,
        item_btwn = item_btwn,
        both_win = both_win))

    # Test whether each subject is only in one of the levels of the between-subject/between-items conditions
    for s in nlevels(subj_n, "S"), cond in [:cond1, :cond2]
        @test length(unique(data[isequal.(data.subj, s), cond])) == 1
    end

    # Test whether each item is only in one of the conditions
    for i in nlevels(item_n, "I"), cond in [:cond1, :cond2]
        @test length(unique(data[isequal.(data.item, i), cond])) == 1
    end

end;

@testset "simdat_crossed test all combinations" begin
    both_btwn = Dict(:cond1 => ["A", "B"],
                     :cond2 => ["C", "D", "E"])
    subj_btwn = merge(Dict(:age => ["O", "Y"]), both_btwn)                    
    item_btwn = merge(Dict(:pet => ["cat", "dog"]), both_btwn)
    both_win = Dict(:time => ["morning", "evening"])


    subj_n = 12
    item_n = 12
    data = DataFrame(simdat_crossed(subj_n, item_n,
        subj_btwn=subj_btwn,
        item_btwn=item_btwn,
        both_win=both_win))
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test all(==("A"), s2.cond1)
    @test all(==("Y"), s2.age)
    @test length(unique(s2.pet)) == 2 # from item between
    
    #-----
    # no subject between
    data = DataFrame(simdat_crossed(subj_n, item_n,
        item_btwn=item_btwn,
        both_win=both_win))

    @test nrow(data) == 288 # many more rows because many more effects are within-subject
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test length(unique(string.(s2.cond1) .* string.(s2.cond2))) == 6 # test all 6 combinations are there
    @test length(unique(s2.cond1)) == 2 # now we have both conditions within subject

    i2 = subset(data, :item => ByRow(==("I02")))
    @test all(==("A"), i2.cond1) # but only one because cond1 is between-items here


    #-----
    # no item between
    data = DataFrame(simdat_crossed(subj_n, item_n,
        subj_btwn=subj_btwn,
        both_win=both_win))
    @test nrow(data) == 288 # many more rows because many more effects are within-subject
    i2 = subset(data, :item => ByRow(==("I02")))
    @test length(unique(string.(i2.cond1) .* string.(i2.cond2))) == 6 # test all 6 combinations are there
    @test length(unique(i2.cond1)) == 2 # now we have both conditions within item
    
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test all(s2.cond1 .== "A") # but only one because cond1 is between-subject here

    #---------
    data = DataFrame(simdat_crossed(subj_n, item_n; both_win))
    @test nrow(data) == subj_n * item_n * 2
    i2 = subset(data, :item => ByRow(==("I02")))
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test length(unique(i2.subj)) == 12
    @test length(unique(s2.item)) == 12


    #------

    data = DataFrame(simdat_crossed(subj_n, item_n))
    @test nrow(data) == subj_n * item_n


    #-----
    data = DataFrame(simdat_crossed(subj_n, item_n; subj_btwn, item_btwn))
    @test nrow(data) == 24
    i2 = subset(data, :item => ByRow(==("I02")))
    @test nrow(i2) == 2
    @test length(unique(i2.subj)) == 2 # because age is within item, but between subjects, we have two subjects here
    @test length(unique(i2.age)) == 2 # because age is within item

    #----
    data = DataFrame(simdat_crossed(subj_n, item_n,
        subj_btwn=subj_btwn))
    @test nrow(data) == 12*12
    i2 = subset(data, :item => ByRow(==("I02")))
    @test length(unique(i2.subj)) == 12
    @test length(unique(i2.cond1 .* i2.cond2 .* i2.age)) == 12 # everything is within item
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test length(unique(i2.subj)) == 12
    @test length(unique(s2.cond1 .* s2.cond2 .* s2.age)) == 1 # everything is within item
    
    #--- same as previous but for item
    data = DataFrame(simdat_crossed(subj_n, item_n,
        item_btwn=item_btwn))
    @test nrow(data) == 12 * 12
    i2 = subset(data, :item => ByRow(==("I02")))
    @test length(unique(i2.subj)) == 12
    @test length(unique(i2.cond1 .* i2.cond2 .* i2.pet)) == 1 # everything is within item
    s2 = subset(data, :subj => ByRow(==("S02")))
    @test length(unique(i2.subj)) == 12
    @test length(unique(s2.cond1 .* s2.cond2 .* s2.pet)) == 12 # everything is within item


end
