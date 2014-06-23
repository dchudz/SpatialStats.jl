using SpatialStats; s = SpatialStats
using Distributions
using Distance

ntrain=50
xtrain = rand(ntrain)
xtest = [0:.01:1]
cov = s.SquaredExponentialCov(length=.2, sigma_spatial = 1, sigma_nonspatial=0.2)
C = s.selfcov(cov, xtrain)


y = rand(MvNormal(C))
yhat = (s.crosscov(cov,xtest,xtrain)/C)*y

using Gadfly
using DataFrames

traindf = DataFrame(x=xtrain, y=y, usage="train")
testdf = DataFrame(x=xtest, y=yhat, usage="test")

p = plot(rbind(traindf, testdf), x=:x, y=:y, color=:usage)
