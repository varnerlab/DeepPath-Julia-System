# include -
include("Include.jl")

# set the path, and load -
path_to_file = joinpath(_PATH_TO_DATA, "KEGGCpdEntry_Smiles.csv");
df = loaddataframe(path_to_file);
compound_dictionary = build(MyKeggChemicalCompoundDictionary(),df)

# get rule: C05378 <=> C00111 + C00118
reactants = [
    compound_dictionary["cpd:C00668"]
];

products = [
    compound_dictionary["cpd:C05345"]
]

(R, P, Δ) = reactionrule(reactants, products);