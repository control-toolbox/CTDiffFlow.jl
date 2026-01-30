using Pkg
Pkg.activate(".")
#Pkg.add("SciMLSensitivity")
#Pkg.add("BenchmarkTools")
#Pkg.add("Enzyme")
#Pkg.add("Mooncake")
#println(pwd())
#println(Pkg.status())
using Markdown
using LinearAlgebra
using Test
using DifferentiationInterface

using ForwardDiff: ForwardDiff
using Enzyme: Enzyme
using Mooncake: Mooncake
using Zygote: Zygote

using OrdinaryDiffEq
#using SciMLSensitivity

include("./fun_examples.jl")
include("../src/CTDiffFlow.jl")
using .CTDiffFlow


tol_error = 2*eps()

# Vectors
println("Automatic differentiation")
println("--------------------------")
# Problems with Enzyme
#Backend = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
Backends = (AutoForwardDiff(), AutoMooncake(), AutoZygote())
#Backends = (AutoForwardDiff(),)
reltol = 1.e-2;
abstol = 1.e-4
tol_error = 10*max(reltol,abstol)
λ = [1.0, 2]
t0 = 0.0;
tf = 1.0;
x0 = [λ[2], 1.0, 1]; 
ivp = ODEProblem(fun_lin, x0, (t0,tf), λ)
algo = Tsit5()
sol = solve(ivp, alg=algo; reltol = reltol, abstol = abstol)
println("Times for the initial flow = ", sol.t)


println("Linear system, with respect to the initial condition")
    for backend in Backends
       println("backend = ", backend)

      # Diff auto
      # Derivative with respect to x0
      λ = [1.0, 2]
      t0 = 0.0;
      tf = 1.0;
      x0 = [λ[2], 1.0, 1]; 
      ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun_lin, t0, x0, tf, λ; backend = backend)
      println("∂x0_flow = ", ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol, print_times=true))
      println("sol_∂xO_flow(tf,λ) = ", sol_∂xO_flow(tf,λ))
      println("ccc", sol_∂xO_flow(tf,λ)-∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
      
      println(@test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))
    end
