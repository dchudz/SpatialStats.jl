using SpatialStats; s = SpatialStats
using Distributions
using Distance
using DataFrames

ntrain=50
xtrain = [[rand(ntrain)] [rand(ntrain)]]

C = s.cov_exp(xtrain, xtrain) + .1*s.diagcov(ntrain)
y = rand(MvNormal(C))

s.variogram(xtrain,y)