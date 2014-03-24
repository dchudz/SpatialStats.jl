using SpatialStats; s = SpatialStats
using Distributions
using Distance
using DataFrames

ntrain=50
xtrain = [[rand(ntrain)] [rand(ntrain)]]
C = s.selfcov(s.SquaredExponentialCov(1,1,.5), xtrain)

y = rand(MvNormal(C))

s.variogram(xtrain,y)