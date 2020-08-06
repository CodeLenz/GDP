using LinearAlgebra
using SparseArrays
using StaticArrays
using ProgressMeter
using WallE

include("Kele.jl")
include("analysis.jl")
include("global.jl")
include("mesh.jl")
include("gmsh.jl")
include("vizinhanca.jl")
include("compliance.jl")
include("filtro_ng.jl")
include("LA_flex.jl")



function Min_Compliance(nx::Int64,ny::Int64,filtro::String="Projecao",
                        Raio::Float64=0.0,vf::Float64=0.4,
                        simp::Float64=3.0,η::Float64=0.5,
                        bolha::Bool=true,
                        niter::Int64=5,
                        arquivo_saida::String="saida.pos")


      # Vamos nos certificar que a opção de filtro é válida
      filtros = Dict("I"=>0, "Projecao"=>1)
      @assert filtro in keys(filtros) "Opção de filtro não é válida"

      # Gera a nossa amiga viga engastada
      nn,coord,ne,connect,nebc,ebc,nnbc,nbc = Cantilever_Beam_Down(nx,ny,100.0)

      # Abre um arquivo do gmsh para escrevermos os resultados
      Inicializa_Malha_Gmsh(arquivo_saida,nn,ne,connect,coord)

      # Vamos gerar um raio
      if Raio==0.0
         Raio = (coord[2,1]-coord[1,1])*1.5
         #println("Calculando o raio automaticamente ",Raio)
      end

      # Monta o operador linear do filtro. Default para Identidade
      MAP = Diagonal(ones(ne))
      if filtro=="Projecao"
         MAP = Mapeamento_x_ρ(Raio,ne,coord,connect)
      end

      # Dados gerais do problema
      E   = 1.0E5
      v   = 0.3
      esp = 0.1

      # Gera um elemento for all. Recupera tanto a matriz de rigidez quanto
      # a matriz de localização da bolha, mesmo se esta não for utilizada
      Ke, Ae = Kele(1,coord,connect,E,v,esp,bolha)

      # Lista de elementos fixos
      lista_fixos = [] #Int.(collect(nx*(ny/2)+1:nx*(ny/2)+20))
      valor_fixos = [] #ones(length(lista_fixos))

      # Variáveis de projeto
      # x = 1.0*ones(ne)
      x = vf*ones(ne)

      # Faz uma cópia das densidades relativasprojetadas
      # só para podermos devolver no final
      ρ_proj = copy(x)

      # Restrições laterais
      ci = zeros(ne)
      cs = ones(ne)

      # Garante que os fixos se mantenham
      ci[lista_fixos] .= valor_fixos
      cs[lista_fixos] .= valor_fixos
      x[lista_fixos] .= valor_fixos

      # Vetores para acompanhar o objetivo e a restrição
      #VF = zeros(niter)
      #VV = zeros(niter)


      # Loop de otimização começa aqui
      niter_penal = 10

      # Se não tivermos filtro, então não precisamos
      # de iterações de continuação (β)
      if filtro=="I"
         niter_penal = 1
      end

      β = 1.0
      β_Δ = 3.0

      # Penalização inicial
      c = 100.0
      mu = zeros(1)


      for iter_proj=1:niter_penal


         # Gera os alias para as chamadas do Driver_Flex_Vol
         f(x) = Driver_Flex_Vol(x,MAP,ne,β,η,
           nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,simp,
           mu,c,vf, lista_fixos, valor_fixos,arquivo_saida,"LA")

         df(x) = Driver_Flex_Vol(x,MAP,ne,β,η,
          nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,simp,
          mu,c,vf,lista_fixos, valor_fixos,arquivo_saida,"dLA")

         g(x) = Driver_Flex_Vol(x,MAP,ne,β,η,
           nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,simp,
           mu,c,vf,lista_fixos, valor_fixos,arquivo_saida,"g")

         proj(x) = Driver_Flex_Vol(x,MAP,ne,β,η,
           nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,simp,
           mu,c,vf,lista_fixos, valor_fixos,arquivo_saida,"proj")



         # Cria o tipo de entrada para o otimizador
         options = GDP.Init()
         options["NITER"]=2000
         options["SHOW"]=false
         

         # Este loop é o do LA
         # ou seja, β deve ser fixo aqui....
         #@showprogress "Otimizando LA $β"
         for LA=1:niter


           # Chama o otimizador caseiro...
           outputs =  GDP.Solve(df, x, ci, cs, options)

           # Extrai o ponto devolvido pelo otimizador
           x_opt = outputs["RESULT"]

           # Calcula as restrições no final do loop interno
           g_atual = g(x_opt)

           if options["SHOW"]
              println("Violações ", g_atual)
           end   

           # Atualiza a penalização e os multiplicadores
           c = min(c * 1.1,100.0)
           mu = Oper.(mu .+ c*g_atual)

           # Calcula o rho projetado e grava
           ρ_proj = proj(x_opt)

           x = copy(x_opt)

      end # LA

      # Atualiza o fator de projeção e segue o baile
      β =  β*β_Δ

   end #iter_proj

   Adiciona_Vista_Escalar_Gmsh(arquivo_saida,"Densidades projetadas",ne, ρ_proj,0.0)

   # Retorna as variáveis de projeto otimizadas e as densidades relativas
   # associadas a este x.
   return ρ_proj

end

