module SpatialStats
using Distance
using Gadfly
using DataFrames

include("covariance.jl")

function variogram(x,y; pairs_per_group = 50)
	D = pairwise(Euclidean(), x', x')
	ysqdiffs = pairwise(Euclidean(), y', y')
	pairsdf = DataFrame(dist = reshape(D,length(D)), 
		ysqdiffs = reshape(ysqdiffs, length(ysqdiffs)))
	sort!(pairsdf, cols=:dist)
	pairsdf[:group] = floor([1:nrow(pairsdf)] / pairs_per_group)
	maxgroup = maximum(pairsdf[:group])
	pairsdf[pairsdf[:group] .== maxgroup, :group] -= 1
	variodf = by(pairsdf, :group, 
		df -> DataFrame(
			mean_sq_diff = mean(df[:ysqdiffs]), 
			dist = median(df[:dist]),
			npairs = size(df,1)
			)
		)
	plot(variodf, x=:dist, y=:mean_sq_diff)
end

# should argument be type or object?
# not sure it works as a type... should maybe make it an object, to be used as the initialization
function fit(C::SpatialCovarianceStructure, x, y)
end

function predict(c::SpatialCovarianceStructure, xtrain, ytrain, xtest) 
end




end # module SpatialStats
