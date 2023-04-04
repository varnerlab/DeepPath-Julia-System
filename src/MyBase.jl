import Base.-

function -(m1::MyKeggChemicalCompoundModel,m2::MyKeggChemicalCompoundModel)
    return m1.morganfingerprint - m2.morganfingerprint
end