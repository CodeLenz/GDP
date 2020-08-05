
# 
# Tests for Algorithm 3
#

using Test

#
# Teste
#
function TESTE()


    function f(x)
        return (x[1]-2)^2 + (x[2]+5)^2
    end
    
    function df(x)
        #sleep(5.0)
        return [ 2*(x[1]-2) ; 2*(x[2]+5) ] 
    end
    
    # Initial point
    x0 = [10.0 ; 10.0]

    # Side constraints
    ci = [3.0 ; -30.0]
    cs = [] 


     # Call optimizer
    options = GDP.Init()
    options["NITER"] = 100
     
     
    output =  GDP.Solve(df, x0, ci, cs, options) 

    # Recovering solution
    x_opt = output["RESULT"]
    flag_converged = output["CONVERGED"]
    opt_norm = output["NORM"]

 
    @show x_opt, flag_converged, opt_norm, output["ITERS"]


end

function BOOTH()

    function df(x)
        df1 = 2*(2*x[2]+x[1]-7)+4*(x[2]+2*x[1]-5)
        df2 = 4*(2*x[2]+x[1]-7)+2*(x[2]+2*x[1]-5)
        return [df1 ; df2]
    end

    # Ponto inicial
    x0 = -5*rand(2)

    # Restrições laterais
    ci = -Inf*ones(2)
    cs =  Inf*ones(2)

    # Call de optimizer
    options = Init()
    options["NITER"] = 1000
    output  = GDP(df, x0, ci, cs, options) 
    
    # The test
    @test isapprox(output["RESULT"],[1.0 ; 3.0],rtol=1E-2)
    @test output["CONVERGED"]
    @show output["ITERS"]
   

end