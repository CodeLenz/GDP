using Test
using GDP


# Run tests
 
println("#"^80, "\n"," "^30,"  Test inputs\n","#"^80)
@time include("tests_inputs.jl")

println("#"^80, "\n"," "^40,"  Test 1\n","#"^80)

@time include("tests_unconstrained.jl")

println("#"^80, "\n"," "^40,"  Test 2\n","#"^80)

@time include("tests_constrained.jl")

#println("#"^80, "\n"," "^40,"  Test 3 \n","#"^80)

#@time include("tests_topopt.jl")
