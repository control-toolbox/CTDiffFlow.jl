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
using SciMLSensitivity

#include("./fun_examples.jl")
#include("../src/CTDiffFlow.jl")
#using .CTDiffFlow

# Automatic differentiation on the flow
# derivatives with respect to x0
function build_∂x0_flow(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff())

    # Jacobian matrix
    function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real};ode_kwargs...)
        function _flow(x0)
            ivp = ODEProblem(rhs, x0, (t0,tf), λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            sol = solve(ivp, alg=algo; ode_kwargs...)
            return sol.u[end]
        end
        return jacobian(_flow,backend,x0)
     end
    return ∂x0_flow
end

function main()
    tol_error = 2*eps()
    # Backends
    # Problems with Enzyme
    #Backend = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
    Backends = (AutoForwardDiff(), AutoMooncake(), AutoZygote())
    reltol = 1.e-3;
    abstol = 1.e-6
    tol_error = 10*max(reltol,abstol)
    #
    # Initial value problem
    λ = [1.0, 2]
    A(λ) = [λ[1] 0 ; 0 λ[2]]
    fun_lin1(x,λ,t) = A(λ)*x
    t0 = 0.0; tf = 1.0;
    x0 = [1., 2.]
    # jacobien of the tlow
    sol_∂xO_flow(tf,λ) = exp(tf*A(λ))

    # Test of automatic differentiation
    for backend in Backends
      println("backend = ", backend)
      #∂x0_flow = CTDiffFlow.build_∂x0_flow(fun_lin1, t0, x0, tf, λ; backend = backend)
      ∂x0_flow = build_∂x0_flow(fun_lin1, t0, x0, tf, λ; backend = backend)
      println("∂x0_flow = ", ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
      println("sol_∂xO_flow(tf,λ) = ", sol_∂xO_flow(tf,λ))
      println(@test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))
      @test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error)
    end
  end

  main()