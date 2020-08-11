

module GDP


  using LinearAlgebra, ProgressMeter, Dates
 
  export Solve, Init

  global VALID_INPUTS 
  VALID_INPUTS = ["NITER","TOL_NORM","ALPHA_0","FACTOR_Z","MIN_STEP","SHOW"]

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

    # Check input keys
    global VALID_INPUTS
    for chave in keys(inputs)
        if !(chave in VALID_INPUTS)
            println("Invalid key in GDP::inputs $chave")
            println("Valid keys are $(VALID_INPUTS)")
            error("Invalid key")
        end
    end

    #
    # First thing is to extract the input parameters 
    #
    nmax_iter  = inputs["NITER"]
    tol_norm   = inputs["TOL_NORM"]
    α_0    = inputs["ALPHA_0"]
    factor_z   = inputs["FACTOR_Z"]
    min_α  = inputs["MIN_STEP"]
    flag_show  = inputs["SHOW"]
    
    
    # Check if the length of x0, ci and cs are the same
    @assert length(x0)==length(ci)==length(cs) "Check_inputs:: length of ci, cs and x0 must be the same"

    # Check if x0 is inside the bounds
    @assert  sum(ci .<= x0 .<= cs)==length(x0) "Check_inputs:: x0 must be inside the bounds ci and cs" 

    # Check if nmax_iter is positive
    @assert  nmax_iter > 0 "Check_inputs:: NITER must be larger than zero "

    # Check if tol_norm is in (0,1)
    @assert 0.0<tol_norm<1.0 "Check_inputs:: TOL_NORM must be in (0,1)"

    # Check if α_0 is > 0 
    @assert α_0>0 "Check_inputs:: ALPHA_0 must be > 0"

    # Check if factor_z is < 1
    @assert 0<factor_z<1 "Solve::Check_inputs:: FACTOR_Z must be in (0,1)"

    # Check if min_α is smaller than α_0 and if it is positive as well
    @assert 0 < min_α < α_0 "Solve::Check_inputs:: MIN_STEP must be in (0,ALPHA_0)"

    # Return input parameters to the caller
    return nmax_iter,tol_norm, α_0, factor_z, min_α, flag_show

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
  # Modified Line Search (Armijo). 
  #
  function Armijo_Projected(f::Function,x0::Array{Float64}, D::Array{Float64}, ci::Array{Float64}, cs::Array{Float64},  c::Float64=0.1,
                            τ::Float64=0.5, α_ini::Float64=10.0, α_min::Float64=1E-12)

    # We start with one function evaluation
    f0 = f(x0)

    # Local vectors
    xn = copy(x0)
    Δx = copy(x0)

    # Initial step
    α = α_ini

    # Flag (success)
    flag_success = false

    # Main Loop
    while true

        # Candidate point (xn)
        xn .= max.(ci,min.(cs,x0 - α*D))

        # Δx
        @. Δx = xn-x0 

        # Left side
        fn = f(xn)

        # Rigth side
        right = f0 + c*dot(D,Δx)

        # First Wolfe condition
        if fn <= right
            flag_success = true
            break
        end #fn <= right

        # Otherwise, decrease step    
        α = α*τ

        # Check for minimum step
        if α<=α_min
            break
        end

    end #while true


    # return 
    return α, flag_success

end #Armijo_Projected


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
  f::Function         -> Objective function \\  
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
function Solve(f::Function,df::Function, x0::Array{Float64}, ci=[], cs=[], inputs=Dict()) 

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
    niter, tol, α0, factor_z, α_min, flag_show = Check_inputs(df, x0, ci, cs, inputs)

    # counts if the number of iterations with step <= α_min
    cont_α_min = 0

    # k means actual iteration and K means next iteration

    # List of variables (auxiliar)
    lvar = 1:length(x0)

    # Steps
    lk = α0
    lK = α0

    # Thetas    
    θk = Inf
    θK = Inf

    # Copy of x0, for subsequent iterations
    xk = copy(x0)
    x = copy(x0)
    xK = copy(x0)

    # Initial gradient
    Dk = df(xk) 
    DK = copy(Dk)

    ########################################### Eq. 36 in the reference paper ####################################
    for i in LinearIndices(DK)
       if (  isapprox(xK[i],ci[i],atol=α0) && Dk[i]>0) || ( isapprox(xK[i],cs[i],atol=α0) && Dk[i]<0)
          Dk[i] =0
       end      
    end  

    ###################################################  STEPS 2,3 and 4 in ALg. 3 ###################################################

    # First step 
    @. xK = xk - α0*Dk

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
    delta_m = Float64[]
    delta_M = Float64[]
    contador = 1
    Prg = Progress(niter, 1, "Minimizing objective function...")
    tempo = @elapsed  begin
    for k=2:niter
 
        # Sensitivity at this point
        DK .= df(xK) 

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
            if (  isapprox(xK[i],ci[i],atol=α0) && DK[i]>0) || ( isapprox(xK[i],cs[i],atol=α0) && DK[i]<0)
               DK[i] =0
            end      
        end    

        ######################################### STEP 6 in Alg. 3 ################################################### 
        T1 = sqrt(1+θk)*lk
        T2 = norm(xK-xk) / (2*norm(DK-Dk))
        α  = min(T1,T2)

        if isnan(α)
            if flag_show
               println("It seems that we are stuck, since Δx=$(norm(xK-xk)) and ΔD=$(norm(DK-Dk)). Skipping the loop." ) 
            end   
            break
        end

        @assert α > 0.0 "GPD::Solve:: step is <=0 ($α)  $T1  $T2 $θk $lk $(norm(xK-xk))" 
        if isinf(α) 
            println("Using Armijo's Backtracking LS since the gradient is constant in this region")
            α, flag_armijo = Armijo_Projected(f,xK, DK, ci, cs, 0.1, 0.5, 10.0, α_min)
            println("Armijo's step ",α)
            if !flag_armijo
               break
            end
        end

        ######################################### STEP 9 in Alg. 3 (Eq. 38) ##########################################
        Alpha_S = Float64[]
        for i in free_x
            if DK[i]>0 && ci[i]>-Inf
               push!(Alpha_S, (xK[i]-ci[i])/DK[i])
            elseif DK[i]<0 && cs[i]< Inf
               push!(Alpha_S, (xK[i]-cs[i])/DK[i])
            end        
        end

        ######################################### STEPS 8 to 15 in Alg. 3 ############################################ 
        aflag = true
        while (aflag)
            if α in Alpha_S
               println("It hapened ! $α")
               α = α*factor_z
            else
               aflag = false
            end
        end        

        
        # next step
        @. x = xK - α*DK

        #
        # Project x and obtains active and inactive sets of variables  
        # Those sets are used to check for first order conditions
        # Steps 16, 17 and 18 in Alg. 3
        #
        active_r, active_r_ci, active_r_cs, free_x = Project!(x,ci,cs) 

        #
        # If the minimum step allowed by the user is used in at least 5 iterations, we can bail outputs
        #
        if α <= α_min
            cont_α_min += 1
         else 
            cont_α_min = 1
         end
 
         if cont_α_min==5
             if flag_show
                println(" Step is too small during 5 consecutive iterations. Skipping")
             end   
             break
         end
 

        # Step 19 in Alg. 3
        θ = α/lk

        # Offsets
        xk .= xK
        xK .= x

        lk = lK
        lK = α

        θk = θK
        θK = θ

        # Compute the change in gradient, to show in the report
        change_DK = norm(DK.-Dk)

        Dk .= DK

         # Fancy report for the mob :)
        if flag_show
            ProgressMeter.next!(Prg; showvalues = [
            (:Iteration,k), 
            (:"Norm at free variables",norm_D), 
            (:"Target norm",tol),
            (:"Current Step",α),
            (:"Current θ",θ),
            (:"Change in design variables",norm(xK.-xk)),
            (:"Change in gradient",change_DK),
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



