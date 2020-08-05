@testset "Test inputs" begin


    # Silly and df(x) - Not used since we are forcing errors

    function df(x)
        one(x)
    end

    # Default inputs
    options = GDP.Init()
    
    #
    # First test - Length of x0,ci and cs must be the same.
    # we are also testing if we can call Solve without input 
    # parameters (Solve must use default parameters)
    #
    x0 = ones(10); ci = zeros(5);  cs = 2*ones(10)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)

    x0 = ones(1); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(5)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)    
 
    #
    # Second test - Check if x0 is inside the bounds
    #
    x0 = -1*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)

    x0 = 5*ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)

    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[5] = 10.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)
   
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    x0[2] = -10.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs)
   

    #
    # Third test - nmax_iter > 0
    #
    options["NITER"]=-1
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["NITER"]=0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
   
    # Restore default parameters
    options = GDP.Init()


    #
    # Fourth test - Check if tol_norm is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["TOL_NORM"]=0.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["TOL_NORM"]=1.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["TOL_NORM"]=10.2
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    
           
    # Restore default parameters
    options = GDP.Init()

    #
    # Sixth test - Check if alpha_0 is > 0.0
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["ALPHA_0"]=-1.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["ALPHA_0"]=0.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    
    # Restore default parameters
    options = GDP.Init()

    #
    # Seventh test - Check if factor_z is in (0,1)
    #
    x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
    options["FACTOR_Z"]=-1.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["FACTOR_Z"]=0.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
    options["FACTOR_Z"]=1.0
    @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)


     # Restore default parameters
     options = GDP.Init()

     #
     # Check if MIN_STEP is in (0,alpha_0)
     #
     x0 = ones(10); ci = zeros(10);  cs = 2*ones(10)
     options["MIN_STEP"]=0.0
     @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
     options["MIN_STEP"]=options["ALPHA_0"]+0.01
     @test_throws AssertionError GDP.Solve(df,x0,ci,cs,options)
 
   
end