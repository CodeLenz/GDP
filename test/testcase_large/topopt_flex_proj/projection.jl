
function Steepest_Descent_Proj(f::Function, df::Function, x0::Array{Float64},
                               ci::Array{Float64}, cs::Array{Float64},
                               flag_show::Bool=false,
                               niter=2000,  tol_norm=1E-6,
                               passo_inicial=5.0, fator_corte = 0.5,
                               passo_minimo = 1E-10)

  # Number of elements
  ne = length(x0)

  # Copy x0 to x1 (internal variable) in order to improve speed
  x1 = copy(x0)

  # Testa para ver se as dimensões das restrições laterais estão OK
  @assert length(ci)==length(x0) "Steepest_Descent_Proj:: size(ci)!=size(x0)"
  @assert length(cs)==length(x0) "Steepest_Descent_Proj:: size(cs)!=size(x0)"

  # Assert that any ci .< cs
  for i in LinearIndices(x0)
          @assert ci[i]<=cs[i] "Steepest_Descent_Proj:: ci must be .<= cs "
  end

  # Initialize the variables
  flag_conv = false
  norma = maxintfloat(Float64)
  passo = passo_inicial

  # Reference value
  f0 = f(x0)
  objetivo_inicial = f0

  # First set of constrained elements (block constrainted)
  Ilast  = Int64[]
  Iblock = Int64[]

  # Main loop.
  #prog = ProgressThresh(passo_minimo, "Internal Loop:")
  nblock_sup = 0
  nblock_inf = 0
  contador = 0
  norma_anterior = norma
  d = zero(x0)
  D = zero(x0)

  # used in Iwall
  le = 1:length(x0)

  number_fletcher = 1

  tempo = @elapsed for iter=1:niter

      # Loop do LS com bloqueio
      flag_conv_interna = false

      # Evaluate the sensitivity by calling f!(D,x)
      tempoD = @elapsed D .= df(x0)

      # Test for any box constraint. Modifies D and x0, as well as 
      # return the set of constrained elements
      Ilast = copy(Iblock)
      Iblock, nblock_inf, nblock_sup = Select_Sets!(D,x0,ci,cs)

      # Norm
      norma_anterior = norma
      norma = norm(D)
      contador += 1

      # Set flag_conv to true
      # and breaks if norm<tol_norm
      if norma<=tol_norm
         flag_conv = true
         break
      end

      # Steepest Descent or conjugate gradient (Fletcher and Reeves)
      # Traditional steepest is chosen for the first iteration OR
      # when there was a change in the set of blocked elements
      # OR when we already have used Fletcher for a large number of 
      # iterations.
      method="Steepest"
      if iter==1 || Iblock != Ilast || mod(number_fletcher,ne)==0
         d .= -D/norma
         number_fletcher = 1
      else
         method="Fletcher"
         number_fletcher += 1
         beta = (norma/norma_anterior)^2
         d .= -D .+ beta*d
         d .= d/norm(d)
      end

      
       # Testa a direção de minimização
       produto_busca = dot(D,d)
       if produto_busca > 0.0
       #   d .= -D/norma
       #   number_fletcher = 1
       #   method="Steepest"
           flag_conv_interna = true
           break

       end

       flag_conv_interna = false

       # Passo no começo da busca
       passo0 = passo

        tempoP = @elapsed while passo > passo_minimo

           # Candidate point and check for bounds.  Also returns the 
           # set of blocked variables
           x1 .= x0 .+ passo*d
           Iwall, flag_direction = Wall!(x1,ci,cs,D,x0)
         
           if !flag_direction
             flag_conv_interna = false
             break
           end


           # We can also evaluate the norm without considering
           # the blocked directions
           free_x = filter(x->!(x in Iwall),le)
           norma = norm(D[free_x])
           if norma<tol_norm
              println("Saindo por norma ", norma)
              flag_conv_interna = true
              flag_conv = true
              break
           end


           # Verifica o valor de f neste novo pto
           f1 = f(x1)
           
           # Condição de descida bem simples. Se melhorou,
           # então aceitamos o ponto. Estou fazendo isto devido
           # ao fato de estarmos projetando as variáveis a cada passo
           # o que muda efetivamente a direção de busca. Conforme já
           # provado (dissertação do Gustavo) estas direções ainda são
           # de minimização.
           if f1<f0

             # Accept the point
             #println("--> ",passo, " ", f1)
             x0 .= x1
             f0 = f1

             # Aumenta o passo ligeiramente, mas só se for menor
             # do que o passo inicial
             if passo0 < passo_inicial
                passo = passo0*2.0
             end

             # Aceita o ponto e sai do loop de LS,
             # atualizando o valor de f0 (referência)
             # e do ponto atual
             #flag_conv_interna = true
             #x0 .= x1
             #f0 = f1
             flag_conv_interna = true
             break

          else

             # Não minimizamos o valor da função. Diminui o passo
             passo = passo*fator_corte

          end #if f1<f0

       end #while

       print("Iteração interna $(iter) | norma $(norma) | metodo $(method)::$(number_fletcher)  | passo $(passo) | $(produto_busca)                                \r")
       #flush(STDOUT)

     # Atualiza a barra de progresso
     #ProgressMeter.update!(prog, passo)

     # Se chegamos aqui, então devemos testar pela
     # convergência da norma.
     if flag_conv
        # Temos uma convergência por norma do gradiente. Podemos sair
        break
     end

    # Se flag_conv_interna for falso, então o LS não teve sucesso
    if !flag_conv_interna 
       if flag_show
           println("\nLS_Proj:: Não tivemos convergência interna - Passo muito pequeno $passo")
       end
       flag_conv = false
       break
   end

  end #iter


  # Tivemos uma solução
  if flag_show
      println("\n********************************************************")
      println("Final do Steepest com projeção")
      println("Objetivo inicial   : ", objetivo_inicial)
      println("Objetivo final     : ", f0)
      if objetivo_inicial!=0.0 && f0!=0.0
         println("% de min.          : ", 100*(f0-objetivo_inicial)/objetivo_inicial )
      end
      println("Bloqueios          : ",nblock_inf," ",nblock_sup)
      println("Numero de iterações: ", contador , " de ",niter)
      println("Tempo total [min]  : ", tempo/60.0)
      println("********************************************************")
  end

  # Returns the solution, convergence flag and last norm
  return x0, flag_conv, norma

end


#
# Selects the set containing the elements that have been blocked
# Returns the set, and the number of elements being blocked (down and up)
# AND MODIFIES the gradient and the design variables
#
function Select_Sets!(D::Array{Float64},x::Array{Float64},
                      ci::Array{Float64},cs::Array{Float64})

      # Empty set
      Iblock = Int64[]; sizehint!(Iblock,length(x))

      # Test for any box constraint
      nblock_sup = 0
      nblock_inf = 0
      @inbounds for i in LinearIndices(D)
          if D[i]>=0.0 &&  x[i]<=ci[i]
             x[i] = ci[i]
             D[i] = 0.0
             nblock_inf += 1
             push!(Iblock,i)
          elseif  D[i]<=0.0 && x[i]>=cs[i]
             D[i] = 0.0
             x[i] = cs[i]
             nblock_sup += 1
             push!(Iblock,i)
          end
      end

      return Iblock, nblock_inf, nblock_sup
    end


    #
    # Applies the "wall" x .= max.(ci,min.(x,cs)) --> Modifies x
    #
    function Wall!(x::Array{Float64},ci::Array{Float64},cs::Array{Float64},
                   D::Array{Float64},x0::Array{Float64})

          # List of Blocked variables
          Iblock_x = Int64[]

          for i in LinearIndices(x)
            if x[i]<=ci[i] 
               x[i] = ci[i]
               push!(Iblock_x,i)
            elseif  x[i]>=cs[i] 
               x[i] = cs[i]
               push!(Iblock_x,i)
            end
          end  

          # Test if the projected direction is a minimizer 
          Δ = x-x0

          direction = dot(Δ,D)
          flag_direction = true
          if direction > 0.0
             flag_direction = false
          end
          #@assert direction < 0.0 "Wall!::projected direction is not a minimizer direction $direction"

          return Iblock_x, flag_direction

    end

