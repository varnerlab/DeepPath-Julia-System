# setup paths 
const _ROOT = pwd();
const _PATH_TO_SRC = joinpath(_ROOT, "src")
const _PATH_TO_DATA = joinpath(_ROOT, "data")

# load external packages
using LinearAlgebra
using RDKitMinimalLib
using Pickle
using CSV
using DataFrames
using DataStructures
using Images
using ImageMagick
using ImageIO
using Plots
using Colors
using DelimitedFiles

# load my codes -
include(joinpath(_PATH_TO_SRC, "Types.jl"))
include(joinpath(_PATH_TO_SRC, "Files.jl"))
include(joinpath(_PATH_TO_SRC, "Factory.jl"))