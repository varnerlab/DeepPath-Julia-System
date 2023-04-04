# include -
include("Include.jl")

# set the path, and load -
path_to_file = joinpath(_PATH_TO_DATA, "KEGGCpdEntry_Smiles.csv");
df = loaddataframe(path_to_file);
compound_dictionary = build(MyKeggChemicalCompoundDictionary(),df)

# get Fp for C00668 -
R = compound_dictionary["cpd:C00668"]
P = compound_dictionary["cpd:C05345"]
Δ = P - R