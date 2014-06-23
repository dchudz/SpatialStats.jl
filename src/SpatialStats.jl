module SpatialStats
using Distance
using Gadfly
using DataFrames
using Optim
using Distributions

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

type FittedSpatialModel
	x_train::Matrix{Float64}
	y::Vector{Float64}
	cov_structure::SpatialCovarianceStructure
end

## Fit and predict methods

function fit(x::Array{Float64}, y::Vector{Float64}, cov_structure::SpatialCovarianceStructure)
	
	function neg_log_lik(param)
    	return -logpdf(MvNormal(selfcov(cov_structure, x)), y)
	end

	function grad!(param, g)
		K = selfcov(cov_structure, x)
		invK = inv(K)
		for i in 1:length(param)
			partialK = partial_selfcov_by_param(cov_structure, x, i)
			# should make this more efficient by only computing the entries we need in the trace
			firstterm = (y'*invK)*partialK*(invK*y) 
			secondterm = trace(invK*partialK)/2
			g[i] = firstterm[1] - secondterm
		end
	end

	function neg_log_lik_with_gradient(g, param::Vector)
		  if !(g === nothing)
		  	grad!(param, g)
		  end
		  return neg_log_lik(param)
	end
	
	param, fval, fcount, converged = fminbox(neg_log_lik_with_gradient, cov_structure.param, lower_constraint(cov_structure), upper_constraint(cov_structure))

end

function predict(model::FittedSpatialModel, xtrain, ytrain) 
end


end # module SpatialStats
