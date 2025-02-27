module GeoArrayOps
using StatsBase: skewness
using OffsetArrays: OffsetMatrix
using ImageFiltering: mapwindow, sobel, imgradients
using Distances: Euclidean, euclidean, evaluate
using PaddedViews: PaddedView
using FillArrays: fill
using StaticArrays: @SMatrix, @MMatrix, SMatrix, MMatrix, MVector
using DataStructures: Deque
using Statistics: median, mean
using ImageCore: scaleminmax, Gray
using LocalFilters

include("utils.jl")
include("pmf.jl")
include("smf.jl")
include("psf.jl")
include("plot.jl")
include("spread.jl")
include("terrain.jl")
include("skew.jl")

export ZevenbergenThorne, Horn, MDG
export pmf, smf, psf
export pssm
export pitremoval
export spread, spread2
export roughness, TRI, TPI, slope, aspect, curvature, hillshade, multihillshade
export skb, skbr

end # module
