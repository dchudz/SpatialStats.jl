abstract SpatialCovarianceStructure 

## SquaredExponentialCov
## 

type SquaredExponentialCov <: SpatialCovarianceStructure
	param::Vector
end

function SquaredExponentialCov(;sigma_spatial=1, length=1, sigma_nonspatial=0)
	SquaredExponentialCov([sigma_spatial, length, sigma_nonspatial])
end

function crosscov(c::SquaredExponentialCov, x1, x2)
	sigma_spatial = c.param[1]
	length = c.param[2]
	D = pairwise(Euclidean(), x1', x2')
	C = sigma_spatial^2 * exp(-(D./length).^2)
	C
end

function selfcov(c::SquaredExponentialCov, x)
	sigma_nonspatial = c.param[3]
	crosscov(c, x, x) + diagm(repeat([sigma_nonspatial^2], outer=[size(x)[1]]))
end

lower_constraint(c::SquaredExponentialCov) = [0., 0., 0.]
upper_constraint(c::SquaredExponentialCov) = [Inf, Inf, Inf]

function partial_selfcov_by_param(c::SquaredExponentialCov, x, param_index::Int)
	sigma_spatial = c.param[1]
	length = c.param[2]
	sigma_nonspatial = c.param[3]
	K = selfcov(c, x)
	D = pairwise(Euclidean(), x', x')
	if param_index == 1
		return (2/sigma_spatial)*K
	end
	if param_index == 2
		return (length)*(D.^2).*K
	end
	if param_index == 3
		return 	diagm(repmat([2sigma_nonspatial], size(x,1)))
	end
	error("Asked to differentiate covariance matrix wrt a bad parameter index")
end



## GammaExponentialCov
## 

type GammaExponentialCov <: SpatialCovarianceStructure
	param::Vector
end

function GammaExponentialCov(;sigma_spatial=1, length=1, gamma = 1, sigma_nonspatial=0)
	GammaExponentialCov([sigma_spatial, length, gamma, sigma_nonspatial])
end

function crosscov(c::GammaExponentialCov, x1, x2)
	sigma_spatial = c.param[1]
	length = c.param[2]
	gamma = c.param[3]
	D = pairwise(Euclidean(), x1', x2')
	C = sigma_spatial^2 * exp(-(D./length).^gamma)
	C
end

function selfcov(c::GammaExponentialCov, x1)
	sigma_nonspatial = c.param[4]
	crosscov(c, x1, x1) + diagm(repeat([sigma_nonspatial^2], outer=[size(x1)[1]]))
end

lower_constraint(c::SquaredExponentialCov) = [0., 0., 0., 0.]
upper_constraint(c::SquaredExponentialCov) = [Inf, Inf, Inf]
