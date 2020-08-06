
    #
    # Monta a matriz de rigidez do elemento isoparamétrico de 4 nós
    #
    function Kele(ele,coord,connect,E,v,esp,bolha::Bool=true)


        # Monta a relação constitutiva do material do elemento
        C = EPT(E,v)

        coord2 = @MMatrix zeros(2,4)

        for i=1:4
            coord2[1,i]=coord[connect[ele,i],1]
            coord2[2,i]=coord[connect[ele,i],2]
        end

        G = [-1.0/sqrt(3.0) 1.0
              1.0/sqrt(3.0) 1.0]

        # Calculo do valor de K local
        Ke = @MMatrix zeros(12,12)
        for s=1:2
            for r=1:2
                B,DJ = Bmatrix(G[s,1],G[r,1],coord2)
                Ke .= Ke .+ transpose(B)*C*B*DJ*esp
            end #r
        end #s

       # Extrai os blocos da matriz de rigidez
       Kaa = Ke[1:8,1:8]
       Kab = Ke[1:8,9:12]
       Kbb = Ke[9:12,9:12]

      A = [I ; -Kbb\Kab']

       # Matriz condensada e matriz de localização
       ifelse(bolha,Kaa .- Kab*(Kbb\Kab'), Kaa), A

end

#
# Específico para um elemento retangular de altura h e largura w
#
function Bmatrix(s::Float64,r::Float64,coord2)
    w = coord2[1,2]-coord2[1,1]
    h = coord2[2,3]-coord2[2,2]
    DJ = (w*h)/4.0

    B = @SMatrix [(s-1)/(2*w)   0   -(s-1)/(2*w)    0   (s+1)/(2*w) 0   -(s+1)/(2*w)    0	-(4*r)/w	0	0	0;
                    0   (r-1)/(2*h) 0   -(r+1)/(2*h)    0   (r+1)/(2*h) 0   -(r-1)/(2*h)  0	 0	0	-(4*s)/h;
                  (r-1)/(2*h)    (s-1)/(2*w)  -(r+1)/(2*h)    -(s-1)/(2*w)    (r+1)/(2*h)  (s+1)/(2*w)  -(r-1)/(2*h)  -(s+1)/(2*w)   0   -(4*r)/w  -(4*s)/h   0]
    return B,DJ
end


function EPT(E::Float64,v::Float64)
    #Otimizando utilizando StaticArrays
    #(E/(1-v^2))*SMatrix{3,3}(1.0,v,0.0,v,1.0,0.0,0.0,0.0,((1.0-v)/2.0))
    (E/(1-v^2))*[1.0 v 0.0;v 1.0 0.0;0.0 0.0 ((1.0-v)/2.0)]

end



function Dofs_ele(ele::Int64,connect::Array{Int64,2})

    # Recupera os nós deste elemento
    no1,no2,no3,no4 = connect[ele,1], connect[ele,2], connect[ele,3], connect[ele,4]

    # Retorna as posições globais deste elemento
    SVector{8}(2*(no1-1)+1,2*(no1-1)+2,2*(no2-1)+1,2*(no2-1)+2,2*(no3-1)+1,2*(no3-1)+2,2*(no4-1)+1,2*(no4-1)+2)

    
end

function Mele(densidade::Float64,w::Float64,h::Float64,espessura::Float64)


  M =  @SMatrix [1.0/9   0   1.0/18   0   1.0/36   0   1.0/18   0 ;
                   0   1.0/9  0   1.0/18   0   1.0/36   0   1.0/18 ;
                 1.0/18  0   1.0/9  0   1.0/18   0   1.0/36   0 ;
                   0  1.0/18   0   1.0/9  0   1.0/18   0   1.0/36 ;
                 1.0/36  0   1.0/18   0   1.0/9  0   1.0/18   0 ;
                   0  1.0/36   0   1.0/18   0   1.0/9  0   1.0/18 ;
                 1.0/18  0   1.0/36   0   1.0/18   0   1.0/9  0 ;
                   0  1.0/18   0   1.0/36   0   1.0/18   0   1.0/9]
  
   return densidade*w*h*espessura*M

end