# GDP
Bounding Box Optimizer for large problems where the optimal solution lies on the boundary. The algorithm is a modified Steepest Descent projecting infeasible variables to the boundary of the feasible design space (defined by the side constraints ci and cs) and with no line-search.  

Example

```julia
    using GDP

    # Gradient of objective function   
    function df(x)
        df1 = 2.0*(x[1]-1)-400*x[1]*(x[2]-x[1]^2)
        df2 = 200.0*(x[2]-x[1]^2)
        return [df1 ; df2]
    end

    # Initial point
    x0 = [0.0 ; 3.0]

    # Side constraints
    ci = [-Inf ; 0.5]
    cs = [0.8 ; Inf] 

    # Call optimizer
    options = GDP.Init()
    options["NITER"] = 10_000
    output = GDP.Solve(df,x0,ci,cs,options)

    # Recovering solution
    x_opt = output["RESULT"]
    flag_converged = output["CONVERGED"]
    opt_norm = output["NORM"]
    number_iter = output["NITER"]

```

Default input options are

```julia
    "NITER"=>1000
    "TOL_NORM"=>1E-6
    "ALPHA_0"=>1E-10
    "FACTOR_Z"=>0.99
    "MIN_STEP"=>1E-12
    "SHOW"=>true

```
where NITER is the maximum number of iterations
      TOL_NORM is the tolerance to assess convergence
      ALPHA_0 is the initial step length (do no change unless you read the theory)
      FACTOR_Z (do no change unless you read the theory)
      MIN_STEP is the minimum step used to skip the loop (if this value is used for the last 5 iterations)
      SHOW prints information ou, if false, runs silently
