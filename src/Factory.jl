function _parse_reaction_phrase(q::Queue{Char}, tmp::Queue{Char}, delim::Char, species::Array{String,1})

    if (isempty(q) == true)
        if (isempty(tmp) == false)
            push!(species, join(tmp));
        end
    else

        # grab a char -
        next_char = dequeue!(q)
        if (next_char == delim)
            
            # ok, so we have delim -
            if (isempty(tmp) == false)
                push!(species, join(tmp));
            end

            # empty tmp -
            empty!(tmp)
        else

            # we do *NOT* have a delim, so save this char
            enqueue!(tmp, next_char);
        end

        # call me -
        _parse_reaction_phrase(q, tmp, delim, species);
    end
end

function _parse_reaction_phrase(phrase::String; delim='+', 
    prefix::Union{Nothing, String} = nothing)::Dict{Int64,String}

    # initialize -
    compounds = Array{String,1}()
    q = Queue{Char}()
    tmp = Queue{Char}()
    species = Dict{Int64,String}()
    counter = 1;

    # build initial q -
    chars = collect(phrase);
    for c ∈ chars
        enqueue!(q,c);
    end

    # call -
    _parse_reaction_phrase(q,tmp, delim, compounds)

    # build dict -
    # append a prefix if we have one
    for compound ∈ compounds
        if (isnothing(prefix) == true)
            species[counter] = compound
        else
            species[counter] = prefix*":"*compound
        end
        counter += 1
    end

    # return -
    return species
end

function _build_list_of_species(reactions::Dict{Int64,MyChemicalReactionModel})::Array{String,1}

    # initialize -
    list_of_species = Array{String,1}()

    # process each reaction -
    for (_, reaction) ∈ reactions
        
        stoichiometry = reaction.stoichiometry;
        for (species,coeff) ∈ stoichiometry
            if (in(species, list_of_species) == false)
                push!(list_of_species, species)
            end
        end
    end

    # return -
    return list_of_species
end

function _build_list_of_reactions(reactions::Dict{Int64,MyChemicalReactionModel})::Array{String,1}

    # initialize -
    list_of_reactions = Array{String,1}();

    # process each reaction -
    for (_, reaction) ∈ reactions

        rnname = reaction.rnnumber
        if (in(rnname, list_of_reactions) == false)
            push!(list_of_reactions, rnname);
        end
    end

    # return -
    return list_of_reactions
end

function _build_stoichiometric_matrix(species::Array{String,1}, reactions::Dict{Int64,MyChemicalReactionModel})::Array{Float64,2}

    # get number of species, and reactions -
    number_of_species = length(species);
    number_of_reactions = length(reactions);
    S = Array{Float64,2}(undef, number_of_species, number_of_reactions);
    fill!(S,0.0);

    # build the stochiometric matrix -
    for i ∈ 1:number_of_species
        
        # get a metabolite -
        metabolite_id = species[i]
        
        for j ∈ 1:number_of_reactions

            # grab the reaction object, and then metabolites dictionary -
            metabolite_dictionary = reactions[j].stoichiometry
            if (haskey(metabolite_dictionary, metabolite_id) == true)
                S[i,j] = metabolite_dictionary[metabolite_id];
            end
        end
    end

    # return -
    return S
end

"""
    _build_stoichiometry_dictionary(reaction::String) -> Dict{String,Float64}
"""
function _build_stoichiometry_dictionary(reaction::String; delim::Char = '+',
    prefix::Union{Nothing, String} = nothing)::Dict{String,Float64}
    
    # initialize -
    d = Dict{String,Float64}()
    stoichiometric_dictionary = Dict{String,Float64}()

    # recursive descent -
    phrases = split(reaction,"<=>")
    LD = _parse_reaction_phrase(string(phrases[1]), delim=delim, prefix=nothing);
    RD = _parse_reaction_phrase(string(phrases[2]), delim=delim, prefix=nothing);

    # convert to dictionary for the output
    for (_,item) ∈ LD
        
        # default -
        word = lstrip(rstrip(item))

        # check: does this item have an internal space?
        if (contains(word, " ") == true)
            components = split(word," ");
            d[string(components[2])] = -1.0*parse(Float64, string(components[1]));
        else
            d[word] = -1.0
        end
    end

    for (_,item) ∈ RD
        
        # default -
        word = lstrip(rstrip(item))

        # check: does this item have an internal space?
        if (contains(word, " ") == true)
            components = split(word," ");
            d[string(components[2])] = 1.0*parse(Float64, string(components[1]));
        else
            d[word] = 1.0
        end
    end

    # ok, so finally -
    if (isnothing(prefix) == true)
        stoichiometric_dictionary = d;
    else

        # update the keys w/a prefix -
        for (key,value) ∈ d
            new_key = "$(prefix):$(key)"
            stoichiometric_dictionary[new_key] = value;
        end
    end

    # return -
    return stoichiometric_dictionary
end

"""
    build(type::Type{MyChemicalReactionModel},data::NamedTuple) -> MyChemicalReactionModel
"""
function build(type::Type{MyChemicalReactionModel},data::NamedTuple; 
    prefix::Union{Nothing, String} = nothing)::MyChemicalReactionModel
    
    # initialize -
    model = MyChemicalReactionModel();
    
    # add stuff to the model -
    model.ecnumber = data[:ecnumber];
    model.rnnumber = data[:rnnumber]
    model.ename = data[:ename]
    model.reaction = data[:reaction]
    model.reversible = data[:reversible]
    model.stoichiometry = _build_stoichiometry_dictionary(data[:reaction], prefix=prefix)
  
    # return model -
    return model;
end

"""
    build(model::MyKeggChemicalCompoundDictionary, data::DataFrame; 
        namekey="KEGG_Entry", datakey="smiles", nbits::Int64 = 512, radius::Int64 = 2) -> Dict{String,MyKeggChemicalCompoundModel}
"""
function build(model::MyKeggChemicalCompoundDictionary, data::DataFrame; 
    namekey="KEGG_Entry", datakey="smiles", nbits::Int64 = 512, radius::Int64 = 2)::Dict{String,MyKeggChemicalCompoundModel}

    # initialize -
    number_of_compounds = nrow(data)
    dict = Dict{String, MyKeggChemicalCompoundModel}()

    # details for the morgan finger print
    fp_details = Dict{String, Any}("nBits" => nbits, "radius" => radius)

    # main loop -
    for i ∈ 1:number_of_compounds

        # get data from df -
        name = data[i, Symbol(namekey)];
        smiles = data[i, Symbol(datakey)];

        # build my object -
        model = MyKeggChemicalCompoundModel();
        model.id = name;
        model.smiles = smiles;
        
        # build morganfingerprint -
        mol = get_mol(smiles);
        mfp = get_morgan_fp(mol, fp_details);
        model.morganfingerprint = parse.(Int64, collect(mfp))

        # store the model -
        dict[name] = model
    end
    
    # return -
    return dict
end

function build(type::Type{MyStoichiometricMatrixModel}, reactions::Dict{Int64,MyChemicalReactionModel};
    lowerbound::Float64 = 100.0, upperbound::Float64 = 100.0)::MyStoichiometricMatrixModel

    # initialize
    model = MyStoichiometricMatrixModel();

    # First, build list of compounds -
    model.species = _build_list_of_species(reactions);
    model.reactions = _build_list_of_reactions(reactions);
    model.S = _build_stoichiometric_matrix(model.species, reactions);
  
    # return -
    return model
end

function build(type::Type{MyCompoundFingerprintMatrixModel}, data::Dict{String,MyKeggChemicalCompoundModel}, 
    compounds::Array{String,1}; nbits::Int64 = 512, radius::Int64 = 2)::MyCompoundFingerprintMatrixModel

    # initialize -
    fingerprint_model = MyCompoundFingerprintMatrixModel(); # build empty model
    number_of_compounds = length(compounds)
    fingerprint_array = Array{Float64,2}(undef, nbits, number_of_compounds);

    # main -
    for i ∈ 1:number_of_compounds
        
        # compound id -
        compound_id = compounds[i];

        # get compound model -
        model = data[compound_id];
        morganfingerprint = model.morganfingerprint;
        for j ∈ 1:nbits
            fingerprint_array[j,i] = convert(Float64, morganfingerprint[j])
        end
    end

    # add data -
    fingerprint_model.nbits = nbits
    fingerprint_model.radius = radius
    fingerprint_model.species = compounds
    fingerprint_model.F = fingerprint_array

    # return 
    return fingerprint_model;
end
