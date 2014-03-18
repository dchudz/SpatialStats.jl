reload("$(pwd())/src/SpatialStats.jl")
using SpatialStats; s = SpatialStats
using Distributions
using Distance

function mydraw(p)
	imgname = tempname()
	img = PNG(imgname, 8inch, 8inch)
	draw(img, p)
	run(`open $imgname`)
end


ntrain=50
xtrain = rand(ntrain)
xtest = [0:.01:1]

C = s.cov_exp(xtrain, xtrain) #+ 10*s.diagcov(ntrain)
y = rand(MvNormal(C))
yhat = (s.cov_exp(xtest,xtrain)/C)*y



D = pairwise(Euclidean(), xtrain', xtrain')
ysqdiffs = pairwise(Euclidean(), y', y')



using Gadfly
using DataFrames

traindf = DataFrame(x=xtrain, y=y, usage="train")
testdf = DataFrame(x=xtest, y=yhat, usage="test")



print(traindf)
p = plot(rbind(traindf, testdf), 
	x=:x, y=:y, color=:usage)
mydraw(p)