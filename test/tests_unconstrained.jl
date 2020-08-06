#
# Function with no side constraints
#
# https://www.sfu.ca/~ssurjano/optimization.html
#
#
@testset "No side constraints" begin


    # First test, optimal solution = (3.0 , 5.0)
    println("\n\n############\n  Test 1.1\n############")
    function f(x) 
       (x[1]-3)^2 + (x[2]-5)^2
    end

    function df(x)
    	df1 = 2*(x[1]-3)
    	df2 = 2*(x[2]-5)
    	return [df1 ; df2]
    end

    # Initial point
    x0 = 10*rand(2)

    # Side constraints
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)


    # Call optimizer
    options = GDP.Init()
    options["NITER"] = 100
    output = GDP.Solve(df,x0,ci,cs,options)


    # The tests
    @test isapprox(output["RESULT"],[3.0 ; 5.0],rtol=1E-2)
    @test output["CONVERGED"]
    println("\nNumber of iterations ",output["ITERS"])

    options["NESTEROV"] = true
    output_nesterov = GDP.Solve(df,x0,ci,cs,options)

    @test isapprox(output_nesterov["RESULT"],[3.0 ; 5.0],rtol=1E-2)
    @test output_nesterov["CONVERGED"]
    println("\nNumber of iterations using Nesterov ",output_nesterov["ITERS"])


    #println("\n","# Result #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")




   println("\n\n############\n  Test 1.2\n############")
   #  Booth
   # x0 = (-5,-5)
   # min (1,3), f = 0.0

    function f(x) 
       (x[1]+2*x[2]-7)^2 + (2*x[1]+x[2]-5)^2 
    end

    function df(x)
        df1 = 2*(2*x[2]+x[1]-7)+4*(x[2]+2*x[1]-5)
        df2 = 4*(2*x[2]+x[1]-7)+2*(x[2]+2*x[1]-5)
        return [df1 ; df2]
    end

    # Initial point
    x0 = -5*rand(2)

    # Side constraints
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Call de optimizer
    options = GDP.Init()
    options["NITER"] = 1000
    output  = GDP.Solve(df,x0,ci,cs,options)
    
    # The test
    @test isapprox(output["RESULT"],[1.0 ; 3.0],rtol=1E-2)
    @test output["CONVERGED"]
    println("\nNumber of iterations ",output["ITERS"])
   
    options["NESTEROV"] = true
    output_nesterov = GDP.Solve(df,x0,ci,cs,options)
    @test isapprox(output_nesterov["RESULT"],[1.0 ; 3.0],rtol=1E-2)
    @test output_nesterov["CONVERGED"]
    println("\nNumber of iterations using Nesterov ",output_nesterov["ITERS"])

    #println("\n","# Result #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")


   println("\n\n############\n  Test 1.3\n############")
   # Beale, 1958
   # x0 = [1, 1]'
   # xo = [3, 0.5]'
   #f(xo) = 0
   
    function f(x) 
       (1.5-x[1]*(1-x[2]))^2+(2.25-x[1]*(1-x[2]^2))^2+(2.625-x[1]*(1-x[2]^3))^2 
    end

       
    function df(x)
        df1 = 2*(2.625-x[1]*(1-x[2]^3))*(x[2]^3-1)+ 2*(2.25-x[1]*(1-x[2]^2))*(x[2]^2-1)+ 2*(1.5-x[1]*(1-x[2]))*(x[2]-1)
        df2 = 6*x[1]*x[2]^2*(2.625-x[1]*(1-x[2]^3))+4*x[1]*x[2]*(2.25-x[1]*(1-x[2]^2))+2*x[1]*(1.5-x[1]*(1-x[2]))
        return [df1 ; df2]
    end

    # Initial point
    x0 = [1.0 ; 1.0]

    # Side constraints
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Call the optimizer
    options = GDP.Init()
    options["NITER"] = 1000
    output = GDP.Solve(df,x0,ci,cs,options)
   
    # The test
    @test isapprox(output["RESULT"],[3.0 ; 0.5],rtol=1E-2)
    @test output["CONVERGED"]
    println("\nNumber of iterations ",output["ITERS"])
    
    
    options["NESTEROV"] = true
    output_nesterov = GDP.Solve(df,x0,ci,cs,options)

    # The test
    @test isapprox(output_nesterov["RESULT"],[3.0 ; 0.5],rtol=1E-2)
    @test output_nesterov["CONVERGED"]
    println("\nNumber of iterations using Nesterov",output_nesterov["ITERS"])
    

    #println("\n","# Result #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")


   println("\n\n############\n  Test 1.4\n############")
   # Goldstein-Price
   # x0 (-2,2) or (2,-2). .
   # min  (0,-1), f =  3.0

    function f(x) 
        ((x[2]+x[1]+1)^2*(3*x[2]^2+6*x[1]*x[2]-14*x[2]+3*x[1]^2-14*x[1]+19)+1)*((2*x[1]-3*x[2])^2*(27*x[2]^2-36*x[1]*x[2]+48*x[2]+12*x[1]^2-32*x[1]+18)+30)
       
    end

       
    function df(x)
        df1 = ((x[2]+x[1]+1)^2*(3*x[2]^2+6*x[1]*x[2]-14*x[2]+3*x[1]^2-14*x[1]+19)+1)*(4*(2*x[1]-3*x[2])*(27*x[2]^2-36*x[1]*x[2]+48*x[2]+12*x[1]^2-32*x[1]+18)+(-36*x[2]+24*x[1]-32)*(2*x[1]-3*x[2])^2)+(2*(x[2]+x[1]+1)*(3*x[2]^2+6*x[1]*x[2]-14*x[2]+3*x[1]^2-14*x[1]+19)+(x[2]+x[1]+1)^2*(6*x[2]+6*x[1]-14))*((2*x[1]-3*x[2])^2*(27*x[2]^2-36*x[1]*x[2]+48*x[2]+12*x[1]^2-32*x[1]+18)+30)
        df2 = ((x[2]+x[1]+1)^2*(3*x[2]^2+6*x[1]*x[2]-14*x[2]+3*x[1]^2-14*x[1]+19)+1)*((2*x[1]-3*x[2])^2*(54*x[2]-36*x[1]+48)-6*(2*x[1]-3*x[2])*(27*x[2]^2-36*x[1]*x[2]+48*x[2]+12*x[1]^2-32*x[1]+18))+(2*(x[2]+x[1]+1)*(3*x[2]^2+6*x[1]*x[2]-14*x[2]+3*x[1]^2-14*x[1]+19)+(x[2]+x[1]+1)^2*(6*x[2]+6*x[1]-14))*((2*x[1]-3*x[2])^2*(27*x[2]^2-36*x[1]*x[2]+48*x[2]+12*x[1]^2-32*x[1]+18)+30)
        return [df1 ; df2]
    end

    # Initial point
    x0 = [1.0 ; -2.0]

    # Side constraints
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Call the optimizer
    options = GDP.Init()
    options["NITER"] = 1000
    
    output = GDP.Solve(df,x0,ci,cs,options)

    # The test
    @test isapprox(output["RESULT"],[0.0 ; -1.0],rtol=1E-2)
    @test output["CONVERGED"]
    println("\nNumber of iterations ",output["ITERS"])

    
    options["NESTEROV"] = true
    output_nesterov = GDP.Solve(df,x0,ci,cs,options)

    # The test
    @test isapprox(output_nesterov["RESULT"],[0.0 ; -1.0],rtol=1E-2)
    @test output_nesterov["CONVERGED"]
    println("\nNumber of iterations using Nesterov ",output_nesterov["ITERS"])


    #println("\n","# Result #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [x_opt [3.0 ; 5.0]])
    println("\n")



   println("\n\n############\n  Test 1.5\n############")
   # Rosenbrook
   # x0  [0,3]
   # min (1,1), f = 0.0

    function f(x) 
        100*(x[2]-x[1]^2)^2+(x[1]-1)^2
    end

       
    function df(x)
        df1 = 2.0*(x[1]-1)-400*x[1]*(x[2]-x[1]^2)
        df2 = 200.0*(x[2]-x[1]^2)
        return [df1 ; df2]
    end

    # X0
    x0 = [0.0 ; 3.0]

    # Side constraints
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

     # Call the optimizer
    options = GDP.Init()
    options["NITER"] = 100_000
    output = GDP.Solve(df,x0,ci,cs,options)

    # The test
    @test isapprox(output["RESULT"],[1.0 ; 1.0],rtol=1E-2)
    @test output["CONVERGED"]
    println("\nNumber of iterations ",output["ITERS"])
    #@show output
 


    options["NESTEROV"] = true
    output_nesterov = GDP.Solve(df,x0,ci,cs,options)
    
 
    # The test
    @test isapprox(output_nesterov["RESULT"],[1.0 ; 1.0],rtol=1E-2)
    @test output_nesterov["CONVERGED"]
    println("\nNumber of iterations using Nesterov ",output_nesterov["ITERS"])
    #@show output
 

    #println("\n","# Result #")
    #show(IOContext(stdout, :compact => false, :limit => false), "text/plain", [output["RESULT"] [1.0 ; 1.0]])
    println("\n")


end # testset