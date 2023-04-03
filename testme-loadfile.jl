# include -
include("Include.jl")

# setup the path -
path_to_reaction_file = joinpath(_PATH_TO_DATA,"Glycolysis.vff")
d = loadreactionfile(path_to_reaction_file)