

# x -> densidades relativas
# simp -> expoente simp

function Global_Stiffness(ne::Int64,connect::Array{Int64,2},coord::Array{Float64,2},
                          nn::Int64, Ke::Array{Float64,2}, ρ::Array{Float64,1},simp::Float64)


    #println(connect,"Conectividade")
    #println(coord,"coordenadas")
    #Define os vetores I, J e V para o montagem sparse
    # O numero de entradas será de ne*(4*4)
    nent=ne*8*8
    I = zeros(Int64,nent)
    J = zeros(Int64,nent)
    V = zeros(nent)


    # Cria contador para I,J e V
    counter =0

    ρmin = 1E-3

    for ele=1:ne

        # Aplica o SIMP
        Ksimp = (ρmin + (1-ρmin)*ρ[ele])^simp * Ke

        # Obtem os gls do elemento
        dofs = Dofs_ele(ele,connect)

        #Loop para montar I,J e V com a informação para esse elemento
        for i=1:8
            dof_i=dofs[i]
            for j=1:8
                dof_j=dofs[j]
                counter +=1
                I[counter] = dof_i
                J[counter] = dof_j
                V[counter] = Ksimp[i,j]
            end #j
        end #if

    end #ele

    #Montar a matrix global sparse e retorna
    return sparse(I,J,V)

end


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


function Load_Vector(nn::Int64,nnbc::Int64,nbc::Array{Float64,2})

    # Declara o vetor de força
    F = zeros(2*nn)

    #Monta o loop para nbc e guarda informações
    for i=1:nnbc

        # Procura o nó, os graus de liberdade e o valor
        node = nbc[i,1]
        dof = nbc[i,2]
        val = nbc[i,3]

        #Sendo o array ebc é a precisão dupla, todas as operações
        #com o nó e o grau de liberdade também será com dupla precisão.
        #Portanto, para ter o acesso no array posição, temos que converter
        # o valor de Integer

        pos = Int(2*(node-1)+dof)

        #Adiciona o valor na posição
        F[pos] += val
    end

    return F
end
