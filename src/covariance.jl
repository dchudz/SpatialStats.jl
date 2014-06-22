
abstract SpatialCovarianceStructure 

## SquaredExponentialCov
## 

type SquaredExponentialCov <: SpatialCovarianceStructure
	sigma_spatial
	length
	sigma_nonspatial
end

function SquaredExponentialCov(; sigma_spatial=1, length=1, sigma_nonspatial=0)
	SquaredExponentialCov(sigma_spatial, length, sigma_nonspatial)
end

function crosscov(c::SquaredExponentialCov, x1, x2)
	D = pairwise(Euclidean(), x1', x2')
	C = c.sigma_spatial^2 * exp(-(D./c.length).^2)
	C
end

function selfcov(c::SquaredExponentialCov, x1)
	crosscov(c, x1, x1) + diagm(repeat([c.sigma_nonspatial^2], outer=[size(x1)[1]]))
end


## GammaExponentialCov
## 

type GammaExponentialCov <: SpatialCovarianceStructure
	sigma_spatial
	length
	gamma
	sigma_nonspatial
end

function GammaExponentialCov(; sigma_spatial=1, length=1, gamma = 1, sigma_nonspatial=0)
	GammaExponentialCov(sigma_spatial, length, gamma, sigma_nonspatial)
end

function crosscov(c::GammaExponentialCov, x1, x2)
	D = pairwise(Euclidean(), x1', x2')
	print(c.gamma)
	print("\n")
	C = c.sigma_spatial^2 * exp(-(D./c.length).^c.gamma)
	C
end

function selfcov(c::GammaExponentialCov, x1)
	crosscov(c, x1, x1) + diagm(repeat([c.sigma_nonspatial^2], outer=[size(x1)[1]]))
end