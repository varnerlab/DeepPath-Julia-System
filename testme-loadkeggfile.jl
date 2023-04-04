# include -
include("Include.jl")

# setup the path -
path_to_reaction_file = joinpath(_PATH_TO_DATA,"KEGG_Reaction_Data.csv")
reaction_dictionary = loadkeggreactionfile(path_to_reaction_file, prefix="cpd")