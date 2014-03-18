module SpatialStats
using Distance

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


end # module SpatialStats
