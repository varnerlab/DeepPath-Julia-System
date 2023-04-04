function reactionrule(reactants::Array{MyKeggChemicalCompoundModel,1}, 
    products::Array{MyKeggChemicalCompoundModel,1})

    # initialize -
    R = nothing;
    P = nothing;
    
    # get dimensions -
    tmp = reactants[1].morganfingerprint;
    number_of_reactants = length(reactants);
    number_of_products = length(products);
    number_of_bits = length(tmp);

    # compute the number of cols -
    number_of_cols = max(number_of_products, number_of_reactants);
    R = Array{Float64,2}(undef, number_of_bits, number_of_cols)
    P = Array{Float64,2}(undef, number_of_bits, number_of_cols)
    
    # default values are zero -
    fill!(R, 0.0);
    fill!(P, 0.0);

    # fill my reactants array -
    for j ∈ 1:number_of_reactants
        
        # get fingerprint
        mdf = reactants[j].morganfingerprint;
        
        for i ∈ 1:number_of_bits
            R[i,j] = mdf[i] 
        end
    end

    # fill my products array -
    for j ∈ 1:number_of_products
        
        # get fingerprint
        mdf = products[j].morganfingerprint;
        
        # fill -
        for i ∈ 1:number_of_bits
            P[i,j] = mdf[i] 
        end
    end

    # compute the rule -
    Δ = P - R;

    # return -
    return (R, P, Δ)
end