#
# Monta a matriz F que mapeia x para ρ, ou seja,
# filtro de densidades padrão
#
function Mapeamento_x_ρ(R::Float64,ne::Int64,coord,connect)

    # Obtem os vizinhos e pesos para o raio/malha
    vizinhos,pesos = Vizinhos_Pesos(R,ne,coord,connect)

    # Vamos montar uma matriz esparsa de Floats64
    VI = Int64[]
    VJ = Int64[]
    VV = Float64[]

    # Para cada elemento, posiciona os vizinhos nas colunas 
    # Vamos aproveitar e já dividir pelo somatório dos pesos
    for ele=1:ne

        # Somatório dos pesos
        somat = sum(pesos[ele])

        # Varre os vizinhos do elemento
        for par in zip(vizinhos[ele],pesos[ele]/somat) 

            push!(VI,ele)
            push!(VJ,par[1])
            push!(VV,par[2])
             
        end # zip
    end # ele

    return sparse(VI,VJ,VV)
    
end


#
# Operador para a projeção tanh. Este operador
# não é linear (depende de ρ), portanto devemos 
# fazer na forma algoritmica
#
function Mapeamento_ρ_tanh(β::Float64,η::Float64,ρ::Array{Float64,1})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    ρ_proj = Array{Float64}(undef,length(ρ))

    # Loop sobre os elementos, calculando a projeção
    for ele in LinearIndices(ρ)

         ρ_proj[ele] = (tanh(β*(ρ[ele]-η))+tanh(β*η))/(tanh(β*η)+tanh(β*(1-η)))

    end

    return ρ_proj

end

#
# Operador de correção de derivada d ρ_proj / d ρ 
# que será utilizado na regra da cadeia para correção de derivadas
#
function dρ_projdρ(β::Float64,η::Float64,ρ::Array{Float64,1})

    # Primeiro alocamos um vetor com a mesma dimensão 
    # de ρ
    operador = Array{Float64}(undef,length(ρ))

    # Loop sobre os elementos, calculando a projeção
    for ele in LinearIndices(ρ)

         operador[ele] = (β*sech(β*(ρ[ele]-η))^2)/(tanh(β*η)+tanh(β*(1-η)))

    end

    
    return Diagonal(operador)

end

