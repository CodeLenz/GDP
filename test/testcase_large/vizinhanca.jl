
# Rotina para calcular os centróides dos elementros
function Centroides(ne,coord,connect)

   # Define a matriz de centróides
   centroides = Array{Float64}(undef,ne,2)

   # Loop pelos elementos
   for ele=1:ne

       x1,y1 = coord[connect[ele,1],:]
       x3,y3 = coord[connect[ele,3],:]

       centroides[ele,1] = (x1 + x3)/2.0
       centroides[ele,2] = (y1 + y3)/2.0

  end

  return centroides

end


# Rotina para determinar a vizinhança no entorno de um elemento
function Vizinhos_Pesos(R,ne,coord,connect)


    # Contador do número máximo de vizinhos
    num_max = 0

    # Calcula os centróides dos elementos
    centroides = Centroides(ne,coord,connect)

    # Aloca a tabela de vizinhanças
    vizinhos = Array{Int64,1}[]
    pesos = Array{Float64,1}[]

     # Para cada elemento
     for ele=1:ne

         # Pega a posição do centróide do elemento central
         cx,cy = centroides[ele,1], centroides[ele,2]

         # Arrays locais
         vizinhos_ele = Int64[]
         pesos_ele = Float64[]

         # Varre a malha e armazena elementos que estejam dentro do raio
         # de varredura
         for j=1:ne

             # Pega a posição do centróide de j
             cjx,cjy = centroides[j,1], centroides[j,2]

             # Calcula a distância entre os centróides
             d = sqrt( (cjx-cx)^2 + (cjy-cy)^2 )

             # Se for <=R, é vizinho
             if d<=R
                 push!(vizinhos_ele,j)
                 push!(pesos_ele, 1.0 - d/R)
             end #if
         end #j

        # Verifica se o elemento tem vizinhos
        @assert length(vizinhos_ele)>0 "Vizinhança:: elemento $ele não tem vizinhos"

        # Guarda o maior número de vizinhos até o momento
        num_max = max(num_max,length(vizinhos_ele))

        # Armazena nos arrays da malha
        push!(vizinhos,vizinhos_ele)
        push!(pesos,pesos_ele)

     end

     # Informa o número máximo de vizinhos
     #println("Número máximo de vizinhos é ",num_max)

     # retorna somente o número máximo de colunas que foi acessado
     return vizinhos, pesos

end
