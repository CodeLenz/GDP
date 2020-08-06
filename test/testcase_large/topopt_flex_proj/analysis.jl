

function Analysis(ne,nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,ρ,simp)

    #Monta vetor global de forca
    F = Load_Vector(nn,nnbc,nbc)

    # Monta matrix global de rigidez
    K = Global_Stiffness(ne,connect,coord,nn,Ke,ρ,simp)

    # Aplicar condições de contorno essencial
    Apply_Ebc!(nebc,ebc,K,F)

    # Usar Cholesky para resolver o sistema linear
    C = cholesky(Symmetric(K))

    # Soluciona o sistema de equações (equilíbrio)
    U = C\F

    return U,F,C

end
