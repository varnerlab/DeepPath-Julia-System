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
            if (next_char != ' ')
                enqueue!(tmp, next_char);
            end
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


"""
    _build_stoichiometry_dictionary(reaction::String) -> Dict{String,Float64}
"""
function _build_stoichiometry_dictionary(reaction::String; delim::Char = '+')::Dict{String,Float64}
    
    # initialize -
    d = Dict{String,Float64}()

    # recursive descent -
    phrases = split(reaction,"<=>")
    LD = _parse_reaction_phrase(string(phrases[1]), delim=delim, prefix="cpd");
    RD = _parse_reaction_phrase(string(phrases[2]), delim=delim, prefix="cpd");

    # convert to dictionary for the output
    for (_,item) ∈ LD
        d[item] = -1.0;
    end

    for (_,item) ∈ RD
        d[item] = 1.0;
    end

    # return -
    return d
end

"""
    build(type::Type{MyChemicalReactionModel},data::NamedTuple) -> MyChemicalReactionModel
"""
function build(type::Type{MyChemicalReactionModel},data::NamedTuple)::MyChemicalReactionModel
    
    # initialize -
    model = MyChemicalReactionModel();
    
    # add stuff to the model -
    model.ecnumber = data[:ecnumber];
    model.rnnumber = data[:rnnumber]
    model.ename = data[:ename]
    model.reaction = data[:reaction]
    model.reversible = data[:reversible]
    model.stoichiometry = _build_stoichiometry_dictionary(data[:reaction])
  
    # return model -
    return model;
end
