load('DimFormulaSMFVectorValuedLevel1WithCharacter.sage')
load('DimFormulaSMFVectorValuedLevel1WithoutCharacter.sage')
load('DimFormulaSMFScalarValuedLevel1WithCharacter.sage')
load('DimFormulaSMFScalarValuedLevel1WithoutCharacter.sage')

def dim_splitting_smf_degree_2_level_1(j,k,e):
    """
    Put everything together
    """
    k = ZZ(k)
    j = ZZ(j) 
    e = ZZ(e)

    if j == 0 and e == 0: 
      return dim_splitting_SV_All_weight(k)
    elif j == 0 and e == 1:
      return dim_splitting_SV_All_weight_charac(k)
    elif j > 0 and e == 0:
      return dim_splitting_VV_All_weight(k,j)
    elif j > 0 and e == 1:
      return dim_splitting_VV_All_weight_charac(k,j)