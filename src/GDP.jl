

module GDP


  using LinearAlgebra, ProgressMeter, Dates


  export Solve, Init



  # 
  # Generate the dictionary with default values (optional arguments)
  # 
  function Init()

    inputs = Dict()
    push!(inputs,"NITER"=>1000)
    push!(inputs,"TOL_NORM"=>1E-6)
    push!(inputs,"ALPHA_0"=>1E-10)
    push!(inputs,"FACTOR_Z"=>0.99)
    push!(inputs,"MIN_STEP"=>1E-12)
    push!(inputs,"SHOW"=>true)
    
    return inputs

   end

#
# Check if the inputs are consistent
#
function Check_inputs(df::Function, x0::Array{Float64}, ci::Array{Float64}, cs::Array{Float64}, inputs::Dict)

    #
    # First thing is to extract the input parameters 
    #
    nmax_iter  = inputs["NITER"]
    tol_norm   = inputs["TOL_NORM"]
    alpha_0    = inputs["ALPHA_0"]
    factor_z   = inputs["FACTOR_Z"]
    min_alpha  = inputs["MIN_STEP"]
    flag_show  = inputs["SHOW"]
    
    
    # Check if the length of x0, ci and cs are the same
    @assert length(x0)==length(ci)==length(cs) "Check_inputs:: length of ci, cs and x0 must be the same"

    # Check if x0 is inside the bounds
    @assert  sum(ci .<= x0 .<= cs)==length(x0) "Check_inputs:: x0 must be inside the bounds ci and cs" 

    # Check if nmax_iter is positive
    @assert  nmax_iter > 0 "Check_inputs:: NITER must be larger than zero "

    # Check if tol_norm is in (0,1)
    @assert 0.0<tol_norm<1.0 "Check_inputs:: TOL_NORM must be in (0,1)"

    # Check if alpha_0 is > 0 
    @assert alpha_0>0 "Check_inputs:: ALPHA_0 must be > 0"

    # Check if factor_z is < 1
    @assert 0<factor_z<1 "Solve::Check_inputs:: FACTOR_Z must be in (0,1)"

    # Check if min_alpha is smaller than alpha_0 and if it is positive as well
    @assert 0 < min_alpha < alpha_0 "Solve::Check_inputs:: MIN_STEP must be in (0,ALPHA_0)"

    # Return input parameters to the caller
    return nmax_iter,tol_norm, alpha_0, factor_z, min_alpha, flag_show

end


# 
# Generate the dictionary with the outputs
# 
function Outputs(x::Array{Float64},flag_conv::Bool,norm_D::Float64,counter::Int64)

    outputs = Dict()
    push!(outputs,"RESULT"=>x) 
    push!(outputs,"CONVERGED"=>flag_conv)
    push!(outputs,"NORM"=>norm_D)
    push!(outputs,"ITERS"=>counter)

    return outputs

end    


  #
  # Given a point, a search direction and a step
  # return the projected point and the list of
  # constrained directions as well as free directions
  #
  function Project!(x::Array{Float64},ci::Array{Float64},cs::Array{Float64})

    # Lists of active costrained directions
    active_r_ci = Int64[]
    active_r_cs = Int64[]
    active_r    = Int64[]
 
    @inbounds for i in LinearIndices(x)

      # Depending on the search direction, we can test for lower OR upper
      # violations. If violated, store in the arrays
      if x[i]<=ci[i]

         # Keep on the boundary
         x[i] = ci[i]

         # Store 
         push!(active_r_ci,i)
         push!(active_r,i)

     elseif x[i]>=cs[i]

           # Keep on the boundary
           x[i] = cs[i]

           # Store 
           push!(active_r_cs,i)
           push!(active_r,i)
       end

    end #i

   # Free design variables
   fdv = filter(x-> !(x in active_r),1:length(x))

   return active_r, active_r_ci, active_r_cs, fdv

end # Project



#
# Main subroutine. This subroutine follows Algorithm 3 in the reference paper
#
#
#
  """
  GDP.Solve  \\
  \\
  Solve the problem \\
  \\
  Min f(x) \\
  where x ∈ ℜ^n and x ∈ [ci, cs].  \\
  \\
  \\
  The inputs for this function are: \\
  df::Function        -> Gradient of f(x)  -> df(x)->Array{Float64,1} \\
  x0::Array{Float64}  -> Initial point \\
  ci::Array{Float64}  -> Lower side constraints \\
  cs::Array{Float64}  -> Upper side constraints \\
  \\
  Optional inputs are informed as a dictionary with default parameters \\
  \\
  push!(inputs,"NITER"=>1000) \\
  push!(inputs,"TOL_NORM"=>1E-6)  \\
  push!(inputs,"LAMBDA_0"=>1E-10) \\
  push!(inputs,"FACTOR_Z"=>0.99) \\
  push!(inputs,"MIN_STEP"=>1E-12) \\
  push!(inputs,"SHOW"=>true) \\
     
  where \\ 
  
  NITER is (maximum) number of iterations\\
  TOL_NORM is used to asses the first order condition for unbounded variables\\
  ALPHA_0 is the initial step length (see reference manuscript)\\
  FACTOR_Z is used to reduce the step length if it is a singular step\\
  If the step is smaller than MIN_STEP for the last 5 iterations, we terminate the subroutine \\   
  SHOW enables a brief report at the end   \\
  
  Outputs are returned as a dictionary with entries \\
  \\ 
  push!(outputs,"RESULT"=>x)  \\
  push!(outputs,"CONVERGED"=>flag_conv) \\
  push!(outputs,"NORM"=>norm_D) \\
  push!(outputs,"ITERS"=>counter) \\
  \\
  where RESULT is the "optimal" point if CONVERGED is true. NORM is the norm of 
  unbounded variables and ITERS is the effective number of iterations. \\
  \\
  """
function Solve(df::Function, x0::Array{Float64}, ci=[], cs=[], inputs=Dict()) 

    # Size of the problem
    n = length(x0)

    # If ci or cs are empty, we pass them to ±∞
    if length(ci)==0 
       ci = -Inf*ones(n)
    end

    if length(cs)==0
      cs = Inf*ones(n)
    end

    # If Dict inputs is empty, we use default parameters
    if isempty(inputs)
       inputs = Init()
    end


    # Test the inputs for any inconsistency and recover 
    # the input parameters  
    niter, tol, l0, factor_z, alpha_min, flag_show = Check_inputs(df, x0, ci, cs, inputs)

    # counts if the number of iterations with step <= alpha_min
    cont_alpha_min = 0

    # k means actual iteration and K means next iteration

    # List of variables (auxiliar)
    lvar = 1:length(x0)

    # Steps
    lk = l0
    lK = l0

    # Thetas    
    Tk = Inf
    TK = Inf

    # Copy of x0, for subsequent iterations
    xk = copy(x0)

    # Initial gradient
    Dk = df(xk) 

    ###################################################  STEPS 2,3 and 4 in ALg. 3 ###################################################

    # First step 
    xK = xk - l0*Dk

    #
    # Project x and obtains active and inactive sets of variables  
    # Those sets are used to check for first order conditions
    # Steps 2, 3 and 4 in Alg. 3
    #
    active_r, active_r_ci, active_r_cs, free_x = Project!(xK,ci,cs) 

    #################################################################################################################################

    

    ###################################################  STEP 5 in ALg. 3 ###################################################
    flag_converged = false
    norm_D = 0.0
    delta_m = []
    delta_M = []
    contador = 1
    Prg = Progress(niter, 1, "Minimizing objective function...")
    tempo = @elapsed  begin
    for k=2:niter
 
        # Sensitivity at this point
        DK = df(xK) 

        ######################################### Verify for fist order conditions ############################################# 
        # Norm at free positions (not bounded)
        norm_D = norm(DK[free_x])

        # Gradients at x[i]==ci[i] 
        delta_m = DK[active_r_ci]

        # Gradients at x[i]==cs[i]
        delta_M = DK[active_r_cs]

        # Increment counter 
        contador += 1

        # Check for first order conditions 
        if  norm_D<=tol && (all(delta_m .>= 0.0)||isempty(delta_m)) && (all(delta_M .<= 0.0)||isempty(delta_M))
            if flag_show 
                printstyled("First order conditions  $(norm_D) < $tol in $(contador) iterations", color=:green)
            end    
            flag_converged = true
            break
        end

        ########################################### Eq. 36 in the reference paper ####################################
        for i in LinearIndices(DK)
            if (  isapprox(xK[i],ci[i],atol=l0) && DK[i]>0) || ( isapprox(xK[i],cs[i],atol=l0) && DK[i]<0)
               DK[i] =0
            end      
        end    

        ######################################### STEP 6 in Alg. 3 ################################################### 
        T1 = sqrt(1+Tk)*lk
        T2 = norm(xK-xk) / (2*norm(DK-Dk))
        passo  = min(T1,T2)

        if isnan(passo)
            if flag_show
               println("It seems that we are stuck, since Δx=$(norm(xK-xk)) and ΔD=$(norm(DK-Dk)). Skipping the loop." ) 
            end   
            break
        end

        @assert passo > 0.0 "GPD::Solve:: step is <=0 ($passo)  $T1  $T2 $Tk $lk $(norm(xK-xk))" 

        ######################################### STEP 9 in Alg. 3 (Eq. 38) ##########################################
        Alpha_S = Float64[]
        for i in free_x
            if DK[i]>0 && ci[i]>-Inf
               push!(Alpha_S, (xk[i]-ci[i])/DK[i])
            elseif DK[i]<0 && cs[i]< Inf
               push!(Alpha_S, (xk[i]-cs[i])/DK[i])
            end        
        end

        ######################################### STEPS 8 to 15 in Alg. 3 ############################################ 
        aflag = true
        while (aflag)
            if passo in Alpha_S
               println("It hapened ! $passo")
               passo = passo*factor_z
            else
               aflag = false
            end
        end        

        
        # next step
        x = xK - passo*DK

        #
        # Project x and obtains active and inactive sets of variables  
        # Those sets are used to check for first order conditions
        # Steps 16, 17 and 18 in Alg. 3
        #
        active_r, active_r_ci, active_r_cs, free_x = Project!(x,ci,cs) 

        #
        # If the minimum step allowed by the user is used in at least 5 iterations, we can bail outputs
        #
        if passo <= alpha_min
            cont_alpha_min += 1
         else 
            cont_alpha_min = 1
         end
 
         if cont_alpha_min==5
             if flag_show
                println(" Step is too small during 5 consecutive iterations. Skipping")
             end   
             break
         end
 

        # Step 19 in Alg. 3
        T = passo/lk

        # Offsets
        xk .= xK
        xK .= x

        lk = lK
        lK = passo

        Tk = TK
        TK = T

        Dk .= DK

         # Fancy report for the mob :)
        if flag_show
            ProgressMeter.next!(Prg; showvalues = [
            (:Iteration,k), 
            (:"Norm at free variables",norm_D), 
            (:"Target norm",tol),
            (:"Current Step",passo),
            (:"Number of variables in the lower bound",length(active_r_ci)),
            (:"First order for lower bound?",all(delta_m .>= -tol)||isempty(delta_m)),
            (:"Number of variables in the upper bound",length(active_r_cs)),
            (:"First order for upper bound?",all(delta_M .<= tol)||isempty(delta_M))],
            valuecolor = :yellow)
        end

    end #k
    end #timing

     # Final report
    if flag_show
        println("\n********************************************************")
        println("End of the main optimization Loop")
        println("Number of variables    : $(n)")
        println("Free variables         : ", length(free_x))
        println("Blocked variables      : ", length(active_r),": ",  length(active_r_ci)," for lower bound ",length(active_r_cs)," for upper bound")
        println("Number of iterations   : ", contador , " of ",niter)
        println("First order conditions : ", flag_converged, " ", all(delta_m .>= -tol)||isempty(delta_m),
        " ", all(delta_M .<=  tol)||isempty(delta_M))
        println("Norm(free positions)   : ", norm_D," Reference ",tol)
        println("Total time             : ", canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64,tempo)))))
        println("********************************************************")
    end



    # Return the dictionary with the solution
    return Outputs(xK,flag_converged,norm_D,contador)

end


end #module
