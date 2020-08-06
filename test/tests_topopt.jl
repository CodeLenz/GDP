#
# Function with many side constraints 
#
@testset "Topopt" begin


    println("\n\n####################################################################################################################\n  
                 Test 3.1 - Topology Optimization (Compliance with Heaviside projection and Augmented Lagrangian) 
               \n####################################################################################################################")
    
    println("Wait...it takes some time to complete. Number of design variables = $(20*12). Solution at the boundary of S")
    println("Topology can be displayed with: gmsh saida.pos")

    include("testcase_large/topopt_flex_proj/main_LA.jl")
    tempo1= @timed begin xs = Min_Compliance(20,12); end

    # Read a reference solution
    using DelimitedFiles
    xr = readdlm("testcase_large/topopt_flex_proj/ref_x.dat")
    #@test isapprox(xs,xr)

  #=
    println("\n\n####################################################################################################################\n  
                 Test 3.2 - Compare time when using WallE to solve 3.1 (GC) 
               \n####################################################################################################################")
    
    println("Wait...it takes some time to complete. Number of design variables = $(20*12). Solution at the boundary of S")
    println("Topology can be displayed with: gmsh saida_WALLE.pos")
    include("testcase_large/main_LA-WALLE.jl")
    tempo2 = @timed begin xswe = Min_Compliance_WALLE(20,12); end
  =#

end #testset