
#
# Gera uma malha com nx × ny elementos e dimensão Lx × Ly
#
function Mesh(Lx::Float64, nx::Int64, Ly::Float64, ny::Int64)

     # Primeiro geramos os nós de baixo para cima, esquerda para a direita
     nn = (nx+1)*(ny+1)
     coord = zeros(nn,2)

     # Dimensões de cada elemento
     dx = Lx / nx
     dy = Ly / ny


     # Gera a matrix de coordenadas nodais
     x = -dx
     y = -dy
     cont = 0
     for j=1:ny+1
         y += dy
         for i=1:nx+1
             x += dx
             cont += 1
             coord[cont,:] = [x y]
         end #i
         x = -dx
     end #j

     # Gera a matrix de conectividades
     ne = (nx)*(ny)
     connect = zeros(Int64,ne,4)
     no1 = 1
     no2 = 2
     no3 = (nx+1)+2
     no4 = no3-1
     cont = 0
     for j=1:ny
         for i=1:nx
             cont += 1
             connect[cont,:] = [no1 no2 no3 no4]
             no1 += 1; no2 += 1; no3 += 1; no4 += 1
         end #i
         no1 += 1; no2 += 1; no3 += 1; no4 += 1
     end #j


    return nn,coord, ne, connect

end


#
# Cria um problema padrão de otimização topológica
#
function Cantilever_Beam_Down(nx::Int64,ny::Int64, forca::Float64)

    # Gera a malha
    Lx = 8.0
    Ly = 5.0
    nn,coord,ne,connect = Mesh(Lx, nx, Ly, ny)

    # Apoios (engaste na esquerda)
    nebc = 2*(ny+1)
    ebc = zeros(Int64,nebc,2)
    no1 = 1
    cont = 0
    for j=1:ny+1
        for i=1:2
            cont += 1
            ebc[cont,:] = [no1 i]
        end
        no1 += nx+1
    end

    # Informação da força concentrada na ponta inferior
    nnbc = 1
    nbc = [nx+1 2 -forca]


    return nn,coord,ne,connect,nebc,ebc,nnbc,nbc

end
