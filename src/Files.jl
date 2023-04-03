"""
    loadreactionfile(path::String) -> Dict{Int64,MyChemicalReactionModel}
"""
function loadreactionfile(path::String)::Dict{Int64,MyChemicalReactionModel}

     # check: is path legit?
    # in production we would check this path, assume ok for now

    # initialize -
    reactions = Dict{Int64, MyChemicalReactionModel}()
    counter = 0;

    # use example pattern from: https://varnerlab.github.io/CHEME-1800-Computing-Book/unit-1-basics/data-file-io.html#program-read-a-csv-file-refactored
    open(path, "r") do io # open a stream to the file
        for line in eachline(io) # read each line from the stream
            
            # Impl me -
            # line is a line from the file  

            # ecnumber::String
            # rnnumber::String
            # ename::String
            # reaction::String
            # reversible::Bool
            # stoichiometry::Dict{String,Float64}

            # A couple of things to think about: 
            # a) ignore the comments, check out the contains function: https://docs.julialang.org/en/v1/base/strings/#Base.contains
            # b) records are comma delimited. Check out the split functions: https://docs.julialang.org/en/v1/base/strings/#Base.split
            # c) from the data in each reacord, we need to build a MyChemicalReaction object. Each reaction object should be stored in the reactions dict with the name as the key
            if (contains(line,"//") == false && isempty(line) == false)

                fields = split(line, ','); # splits around the ','

                # grab the fields -
                ecnumber = string(fields[1]);
                rnnumber = string(fields[2]);
                ename = string(fields[3]);
                reaction_string = string(fields[4]);
                reversible = parse(Bool, fields[5]);


                # build the reaction model -
                model = build(MyChemicalReactionModel, (
                    ecnumber = ecnumber,  
                    rnnumber = rnnumber,
                    ename = ename,
                    reaction = reaction_string, 
                    reversible = reversible
                ));

                # store -
                reactions[counter] = model;

                # update the counter -
                counter = counter + 1
            end
        end
    end

    # return -
    return reactions;
end