reload("$(pwd())/src/SpatialStats.jl")
using SpatialStats; s = SpatialStats
using Distributions
using Distance
using DataFrames

ntrain=50
xtrain = [[rand(ntrain)] [rand(ntrain)]]

C = s.cov_exp(xtrain, xtrain) + .1*s.diagcov(ntrain)
y = rand(MvNormal(C))

mydraw(s.variogram(xtrain,y))

# mydraw(plot(x=xtrain[:,1], y=xtrain[:,2], color=y))
# mydraw(plot(pairsdf, x=:dist, Geom.histogram))