using IntervalWavelets
using Wavelets
using Base.Test

println("Testing dyadic rationals...")

#=
Translation between indices and values for dyadic rationals
=#

p = 2
S = DaubSupport(-p+1, p)

# Level 0 = integers
supportx = dyadic_rationals( S, 0 )
for x in round(Int, supportx)
	@test supportx[ x2index(x,S) ] == x
end

for idx in 1:length(supportx)
	@test supportx[idx] == index2x(idx,S)
end

# Higher level
R = rand( 1:10 )

supportx = dyadic_rationals( S, R )
for x in supportx
	@test supportx[ x2index(x,S,R) ] == x
end

for idx in 1:length(supportx)
	@test supportx[idx] == index2x(idx,S,R)
end

