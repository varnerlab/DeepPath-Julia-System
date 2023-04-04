abstract type AbstractChemicalReaction end
abstract type AbstractKeggCompoundModel end

struct MyKeggChemicalCompoundDictionary 
end

mutable struct MyCompoundFingerprintMatrixModel

    # data -
    F::Array{Float64,2}
    species::Array{String,1}
    nbits::Int64
    radius::Int64

    # constructor -
    MyCompoundFingerprintMatrixModel() = new()
end

mutable struct MyStoichiometricMatrixModel

    # data -
    S::Array{Float64,2}
    bounds::Array{Float64,2}
    species::Array{String,1}
    reactions::Array{String,1}

    # constructor -
    MyStoichiometricMatrixModel() = new()
end

"""
    MyChemicalReaction <: AbstractKeggReaction

Holds chemical reaction information. See the Test.net file
"""
mutable struct MyChemicalReactionModel <: AbstractChemicalReaction
    
    # data -
    ecnumber::String
    rnnumber::String
    ename::String
    reaction::String
    reversible::Bool
    stoichiometry::Dict{String,Float64}

    # constructor
    MyChemicalReactionModel() = new()
end

"""
MyKeggChemicalCompoundModel
"""
mutable struct MyKeggChemicalCompoundModel <: AbstractKeggCompoundModel
    # data -
    id::String
    smiles::String
    morganfingerprint::Array{Int64,1}

    # constructor
    MyKeggChemicalCompoundModel() = new();
end