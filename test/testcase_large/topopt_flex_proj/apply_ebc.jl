function Apply_Ebc!(nebc::Int64,ebc::Array{Int64,2},K,F)

    #Loop para condição natural de contorno
    for i=1:nebc

        # encontra o no e o grau de liberdade
        no = ebc[i,1]
        dof= ebc[i,2]

        # valida a posicao global
        pos = Int(2*(no-1)+dof)

        # Aplica 0.0 para linha e coluna (equação)
        K[pos,:] .= 0.0
        K[:,pos] .= 0.0
        F[pos] = 0.0

        #Cria diagonal 1.0
        K[pos,pos] = 1.0

    end


end
