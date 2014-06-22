using SpatialStats; s = SpatialStats
using Distributions
using Distance

ntrain=50
xtrain = rand(ntrain)
xtest = [0:.01:1]

cov = s.GammaExponentialCov(gamma = 3, sigma_nonspatial=.9)
#cov = s.SquaredExponentialCov(1,1,.2)
C = s.selfcov(cov, xtrain) #+ 10*s.diagcov(ntrain)
y = rand(MvNormal(C))
yhat = (s.crosscov(cov,xtest,xtrain)/C)*y

using Gadfly
using DataFrames

traindf = DataFrame(x=xtrain, y=y, usage="train")
testdf = DataFrame(x=xtest, y=yhat, usage="test")

p = plot(rbind(traindf, testdf), x=:x, y=:y, color=:usage)

using Distributions
using Optim    

# can return negative numbers, since the optimization is unconstrained and +/- mean the same thing in the cov function

# o = optimize(
# 	pars -> -logpdf(MvNormal(s.selfcov(s.GammaExponentialCov(pars...), xtrain)), y),
# 	[1.,1.,1.,1.]
# 	)

# o = optimize(
# 	pars -> -logpdf(MvNormal(s.selfcov(s.SquaredExponentialCov(pars...), xtrain)), y),
# 	[2.,5.,.3]
# 	)
