# include -
include("Include.jl")

# setup the path -
path_to_reaction_file = joinpath(_PATH_TO_DATA,"Glycolysis.vff")
reaction_dictionary = loadreactionfile(path_to_reaction_file, prefix="cpd")
stm_model = build(MyStoichiometricMatrixModel, reaction_dictionary);
S = stm_model.S;

path_to_file = joinpath(_PATH_TO_DATA, "KEGGCpdEntry_Smiles.csv");
df = loaddataframe(path_to_file);
compound_dictionary = build(MyKeggChemicalCompoundDictionary(),df)
fingerprint_model = build(MyCompoundFingerprintMatrixModel, compound_dictionary, stm_model.species);