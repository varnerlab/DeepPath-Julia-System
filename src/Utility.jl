function find(reactions::Dict{Int64,MyChemicalReactionModel},data::String; 
    key::Symbol = :rnnumber)::Union{Nothing,MyChemicalReactionModel}

    target = nothing;
    for (_,model) ∈ reactions
        
        if (getproperty(model,key) == data)
           target = model;
           break; 
        end
    end

    # return -
    return target
end