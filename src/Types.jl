abstract type AbstractChemicalReaction end

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