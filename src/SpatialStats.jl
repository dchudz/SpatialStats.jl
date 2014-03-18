module SpatialStats
using Distance
using Gadfly
using DataFrames

function cov_exp(x1,x2)
	D = pairwise(Euclidean(), x1', x2')
	C = exp(-(abs(D)))
	C
end

function cov_sq_exp(x1,x2)
	D = pairwise(Euclidean(), x1', x2')
	C = exp(-(100 * D.^2))
	C
end

function diagcov(n)
	diagm(repeat([.01], inner=[n]))
end


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
	p = plot(variodf, x=:dist, y=:mean_sq_diff)
	p
end



end # module SpatialStats
