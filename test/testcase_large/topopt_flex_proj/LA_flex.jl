
# Operador <>
function Oper(a::Float64)
   max(a,0.0)
end

#
#
# Rotina driver para o problema de minimização de Flex
# com restrição de volume
#
#
function Driver_Flex_Vol(x::Array{Float64},MAP,ne,β,η,
                         nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,simp,
                         mu,c,Vlim, lista_fixos, valor_fixos,
                         arquivo_saida::String, opcao::String)


    @assert opcao in ["LA","dLA","g","proj"]

    # O cálculo do equilíbrio e do volume são 
    # comuns a todas as opções
    
    # Mapeamento para as densidades - Filtro
    ρ = MAP*x
     
    # Projeção heaviside
    ρ_proj = Mapeamento_ρ_tanh(β,η,ρ)
   
    if opcao=="proj"
      return ρ_proj
    end

    # Calcula o equilíbrio da estrutura
    U,F,_ = Analysis(ne,nn,nebc,nnbc,connect,coord,ebc,nbc,Ke,ρ_proj,simp)
                 
    # Podemos calcular tanto o objetivo quanto
    # a restrição
    objetivo = Compliance(U,F)
    g_volume = [(Volume(ne,ρ_proj)/Vlim - 1.0)]

    # Dependendo da opção, podemos calcular somente
    # o necessário
    #
    if opcao=="LA"
 
       return objetivo + (c/2.0)*sum( Oper.(mu/c + g_volume).^2 )

    elseif opcao=="g"

    	return g_volume

    elseif opcao=="dLA"

       # Calcula as derivadas w.r.t ρ_proj
       dflex = DCompliance(ne,Ke,coord,connect,U,ρ_proj,simp)
       dvol  = dVolume(ne)

       # Operador de correção da projeção
       R = dρ_projdρ(β,η,ρ) 

       # Corrige as derivadas devido ao filtro
       # Cuidado com o MAP' 
       dflex_corr = R*MAP'*dflex
       dvol_corr  = R*MAP'*dvol

       # Garante os valores fixos
       dflex_corr[lista_fixos] .= 0.0
       dvol_corr[lista_fixos] .= 0.0

       # Monta a derivada do LA
       dLA = dflex_corr .+ Oper.(mu .+ c*g_volume)*(1.0/Vlim).*dvol_corr

       return dLA 
 
    end 	
  

end