function makelevels(nlev, tag='S', inner=1, outer=1)
    levs = string.(tag, lpad.(collect(1:nlev), ndigits(nlev), '0'))
end
