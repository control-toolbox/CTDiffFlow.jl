using Pkg
Pkg.activate(".")
#Pkg.add("ReverseDiff")
#Pkg.add("DataFrames")
#Pkg.add("SciMLSensitivity")
#Pkg.add("BenchmarkTools")
#Pkg.add("Enzyme")
#Pkg.add("Mooncake")
#println(pwd())
#println(Pkg.status())
using DataFrames
using Markdown
using LinearAlgebra
using Test
using DifferentiationInterface

using ForwardDiff: ForwardDiff
#using Enzyme: Enzyme
using Mooncake: Mooncake
#using Zygote: Zygote
#using ReverseDiff: ReverseDiff

using OrdinaryDiffEq


#include("./fun_examples.jl")
include("../src/CTDiffFlow.jl")
using .CTDiffFlow
include("../src/myode43/myode43.jl")


#
#
# Initial value problem
λ = [1.0, 2]
A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
fun_lin(x,λ,t) = A(λ)*x
t0 = 0. ; tf = 1.
tspan = (t0,tf)
x0 = [1., 2., 3]
# jacobien of the flow
sol_∂xO_flow = exp(tf*A(λ))



function test_algo_ode(fun::Function, tspan::Tuple{<:Real,<:Real}, x0::Vector{<:Real}, λ::Vector{<:Real},
         sol_∂xO_flow::Matrix{<:Real},  
         Algorithmes = (("myode43", 7), (RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)), 
         Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote()),
         reltol = 1.e-3, abstol = 1.e-6)
  """
    fun, tspan, x0 and λ : definition of the ode
    Tests of different algorithmes of numerical integration 
    with different backends of automatic differentiation
    input
    -----
    adaptive : Boolean
               true  : variable steps
               false : fixed steps
    Algorithmes : Tuples
                  Algorithmes[i] : Tuples = name of an algorithm, number of caracters for print
    Backends : Tuple
               list of backends of automatic differentiation
    return
    ------
    df_ode : dataframe
  """
  
    tol_error = 2*eps()
    tol_error = 10*max(reltol,abstol)
    t0 = tspan[1]; tf = tspan[2];


    # Test of convergence with different the numerical integration algorithms
    df_algo = DataFrame(VAR_IND=String[], backend=String[], adaptive=Bool[], cv=Bool[], retcode=[])
    for adaptive in (true, false)
    for algorithme in Algorithmes
        algo = algorithme[1]
        name_algo = string(algorithme[1])[1:algorithme[2]]
        # Integration of the IVP
        if algo == "myode43"
        T,X = myode43(fun,x0,λ,(t0,tf),reltol,abstol)
        push!(df_algo, ["IVP", name_algo, adaptive, true, true])
        else 
        ivp = ODEProblem(fun, x0, (t0,tf), λ)
        sol = solve(ivp, alg=algo, reltol = reltol, abstol = abstol)
        push!(df_algo, ["IVP", name_algo, adaptive, true, sol.retcode])
        T = sol.t
        end
        for var_ind in ("IND", "VAR2", "VAR1")
          for backend in Backends
            if var_ind == "IND"
              RelTol = reltol
              AbsTol = abstol
              ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ; var_ind=:ind, backend = backend)
            elseif var_ind =="VAR2"
              RelTol = reltol
              AbsTol = abstol
              ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ; backend = backend)
            else
              my_Inf = prevfloat(typemax(Float64))
              n = length(x0)
              p = n
              RelTol = reltol*ones(n,n+1)
              #RelTol = reltol*ones(n,n) # ==> error of dimension
              AbsTol = abstol*ones(n,n+1)
              ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ; backend = backend)
            end
              try
                  sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, adaptive=adaptive, reltol=RelTol, abstol=AbsTol)
                  push!(df_algo, [var_ind, string(backend), adaptive, true,true])
              catch
                  push!(df_algo, [var_ind, string(backend), adaptive, false,false])
              end
          end
        end
    end 
  end
    return df_algo
  end

df_algo = test_algo_ode(fun_lin,tspan,x0,λ,sol_∂xO_flow)
println(df_algo)
