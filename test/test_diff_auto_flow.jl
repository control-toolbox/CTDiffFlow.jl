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
using Zygote: Zygote
#using ReverseDiff: ReverseDiff

using OrdinaryDiffEq


#include("./fun_examples.jl")
include("../src/CTDiffFlow.jl")
using .CTDiffFlow
include("../src/myode43/myode43.jl")

function main(adaptive)
    tol_error = 2*eps()
    # Backends
    # Problems with Enzyme and Zygote
    #Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
    #Backends = (AutoForwardDiff(), AutoMooncake())
    Backends = (AutoForwardDiff(),)
    #Backends = (AutoMooncake(),)
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

    λ = [1.0, 2]
    B(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
    fun_lin2(x,λ,t) = B(λ)*x
    x0 = [1., 2., 3]
    # jacobien of the tlow
    sol_∂xO_flow2(tf,λ) = exp(tf*B(λ))

    df_sol = DataFrame(VAR_IND=String[], backend=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[], length_times=Real[])
    Sol = []
    
    #algo = RK4()
    #algo = Tsit5()
    algo = "myode43"
    # Integration of the IVP
    if algo == "myode43"
      T,X = myode43(fun_lin2,x0,λ,(t0,tf),reltol,abstol)
      push!(df_sol, ["IVP", "", NaN, NaN, T[2:3], length(T)])
    else
      ivp = ODEProblem(fun_lin2, x0, (t0,tf), λ)
      sol = solve(ivp, alg=Tsit5(), reltol = reltol, abstol = abstol)
      push!(df_sol, ["IVP", "", NaN, NaN, sol.t[2:3], length(sol.t)])
    end
    


    N = 10
    dt = (tf-t0)/N
    # Test of automatic differentiation

      ind = 1
      #for var_ind in ("IND", "VAR2", "VAR1")
      for var_ind in ("IND",)
        for backend in Backends
        if var_ind == "IND"
          RelTol = reltol
          AbsTol = abstol
          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun_lin2, t0, x0, tf, λ; var_ind=:ind, backend = backend)
        elseif var_ind =="VAR2"
          RelTol = reltol
          AbsTol = abstol
          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun_lin2, t0, x0, tf, λ; backend = backend)
        else
          my_Inf = prevfloat(typemax(Float64))
          n = length(x0)
          p = n
          RelTol = reltol*ones(n,n+1)
          #RelTol = reltol*ones(n,n) # ==> error of dimension
          AbsTol = abstol*ones(n,n+1)

          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun_lin2, t0, x0, tf, λ; backend = backend)
        end
        if adaptive
            sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, reltol=RelTol, abstol=AbsTol)
        else
            sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, adaptive=false, dt=dt, reltol=RelTol, abstol=AbsTol)
        end
        push!(Sol,sol)

      #println(sol)
      #println("sol_∂xO_flow2(tf,λ) = ", sol_∂xO_flow2(tf,λ))
        norm_inf = norm(sol-sol_∂xO_flow2(tf,λ),Inf)
        if ind==1
            norm_diff = NaN
        else
            norm_diff = norm(Sol[ind]-Sol[ind-1],Inf)
        end
        ind = ind+1
        push!(df_sol, [var_ind, string(backend), norm_inf, norm_diff, T[2:3], length(T)])
      #println("∂x0_flow = ", ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
      #println("sol_∂xO_flow(tf,λ) = ", sol_∂xO_flow(tf,λ))
      #println(@test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))
       #  @test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error)
      end
    end

    return df_sol,Sol
  end

#df_sol, Sol = main(false)
#println(df_sol)

df_sol, Sol = main(true)
println(df_sol)

using SciMLSensitivity
df_sol, Sol = main(true)
println(df_sol)

