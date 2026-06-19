using Pkg
Pkg.activate(".")
#Pkg.add("Aqua")
using Test
using Aqua
using CTDiffFlow
using LinearAlgebra

using CTFlows
using OrdinaryDiffEq
using ADTypes  # pour ad_backend=ADTypes.AutoForwardDiff()

using DifferentiationInterface
using ForwardDiff: ForwardDiff

include("ode_examples.jl")
#
@testset verbose = true showtiming = true "CTDiffFlow tests" begin
    #for name in (:aqua, :default, :eq_var_rhs)
    for name in (:eq_var_rhs, :end_ind_var, :CTFlows)
    #for name in (:CTFlows,)
        @testset verbose = true "$(name)" begin
            test_name = Symbol(:test_, name)
            include("$(test_name).jl")
            @eval $test_name()
        end
    end
end
println()
