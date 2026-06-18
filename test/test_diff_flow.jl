using Pkg
Pkg.activate(".")
#Pkg.add("DifferentialEquations")
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
#using Mooncake: Mooncake
#using Zygote: Zygote
#using ReverseDiff: ReverseDiff

using OrdinaryDiffEq

using OrdinaryDiffEqLowOrderRK

using OrdinaryDiffEqHighOrderRK

using OrdinaryDiffEqFIRK

include("../src/CTDiffFlow.jl")
using .CTDiffFlow
#include("./fun_examples.jl")

function main(adaptive; 
    internalnorm = :default,
    wrt = :x0,
      #Algorithmes = (Euler(), (RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9))
    Algorithmes = ((Tsit5(), 5),)
    )
    println("wrt = ", wrt)
    println("adaptive = ", adaptive)
    tol_error = 2*eps()
    reltol = 1.e-3;
    abstol = 1.e-6
    tol_error = 10*max(reltol,abstol)
    #
    # Initial value problem
    λ = [1.0, 2]
    A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
    fun_lin1(x,λ,t) = A(λ)*x
    t0 = 0. ; tf = 1.
    x0 = [1., 2., 3]
    # jacobien of the flow
    xf = exp((tf-t0)*A(λ))*x0
    if wrt == :x0
      sol_∂flow = exp((tf-t0)*A(λ))
    elseif wrt == :λ
      sol_∂flow = (tf-t0)*[xf[1] 0 ; 0 xf[2] ; xf[3] -xf[3]]
    elseif wrt == :t0 
      sol_∂flow = [ -xf[1]*λ[1] , -xf[2]*λ[2] , xf[3]*(λ[2]-λ[1])]
    elseif wrt ==:tf
      sol_∂flow = fun_lin1(xf,λ,tf)
    end


    # Test of convergence with different the numerical integration algorithms

   
    df_algo = DataFrame(VAR_IND=String[], backend=String[], adaptive=Bool[], cv=Bool[])
    

    #Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
    Backends = (AutoForwardDiff(),)
      

# Test of automatic differentiation
# ---------------------------------


    df_sol = DataFrame(VAR_IND=String[], backend=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[])
    Sol = []
    ind = 1
    for algorithme in Algorithmes
      algo = algorithme[1]
      name_algo = string(algorithme[1])[1:algorithme[2]]
      # Integration of the IVP
        ivp = ODEProblem(fun_lin1, x0, (t0,tf), λ)
        sol = solve(ivp, alg=algo, reltol = reltol, abstol = abstol)
        push!(df_sol, ["IVP", name_algo, NaN, NaN, sol.t[2:3]])
        T = sol.t

      for var_ind in ("IND", "VAR", "END")
      #for var_ind in ("END",)
        for backend in Backends
          if var_ind == "IND"
            RelTol = reltol
            AbsTol = abstol
            ∂_flow = CTDiffFlow.build_∂flow(fun_lin1, t0, x0, tf, λ; var_ind = :ind, wrt = wrt, internalnorm = internalnorm, backend = backend)
          elseif var_ind == "VAR"
            ∂_flow = CTDiffFlow.build_∂flow(fun_lin1, t0, x0, tf, λ; wrt = wrt, backend = backend)
          elseif var_ind == "END"
            ∂_flow = CTDiffFlow.build_∂flow(fun_lin1, t0, x0, tf, λ; var_ind= :end, wrt = wrt, δh=1.e-10)
          end
          if adaptive
              sol = ∂_flow(t0, x0, tf, λ; alg=algo, adaptive=true, reltol=reltol, abstol=abstol)
          else
            N = 10
            dt = (tf-t0)/N
            sol = ∂_flow(t0, x0, tf, λ; adaptive=adaptive, dt=dt, alg=algo, reltol=reltol, abstol=abstol)
          end
          push!(Sol,sol)
          norm_inf = norm(sol-sol_∂flow,Inf)
          if ind==1
            norm_diff = NaN
          else
            norm_diff = norm(Sol[ind]-Sol[ind-1],Inf)
          end
          ind = ind+1
          push!(df_sol, [var_ind, string(backend), norm_inf, norm_diff, T[2:3]])
        end
      end
    end
    return df_algo, df_sol,Sol
  end

#=
# in the automatic differentiation of the flow there is h'(p) the step derivative, 
# so the diagram doesn't switch
# my_norm2 is the norm of 
# https://github.com/ODINN-SciML/DiffEqSensitivity-Review/blob/main/code/SensitivityForwardAD/example-AD-tolerances.jl
sse2(x::Number) = x^2
sse2(x::ForwardDiff.Dual) = sse2(ForwardDiff.value(x)) + sum(sse2, ForwardDiff.partials(x))
totallength2(x::Number) = 1
function totallength2(x::ForwardDiff.Dual)
  totallength2(ForwardDiff.value(x)) + sum(totallength2, ForwardDiff.partials(x))
end
totallength2(x::AbstractArray) = sum(totallength2, x)
my_norm2 = (u, t) -> sqrt(sum(x -> sse2(x), u) / totallength2(u))
df_algo2, df_sol2, Sol2 = main(true, my_norm2)
#println(df_algo)
println(df_sol2)
=#

#=
# with my_norm the diagram switches 
sse(x::Number) = x^2
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) #+ sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
function totallength(x::ForwardDiff.Dual)
  totallength(ForwardDiff.value(x)) #+ sum(totallength, ForwardDiff.partials(x))
end
totallength(x::AbstractArray) = sum(totallength, x)
my_norm = (u, t) -> sqrt(sum(x -> sse(x), u) / totallength(u))

df_algo, df_sol, Sol = main(true,my_norm)
#println(df_algo)
println(df_sol)
=#
df_algo, df_sol, Sol = main(true,wrt=:tf, Algorithmes = ((Tsit5(), 5),)) #, Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)) )
#println(df_algo)
println(df_sol)
for i in 1:length(Sol)
  println(Sol[i])
end

#=
df_algo, df_sol, Sol = main(false,wrt=:x0, Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)) )
#println(df_algo)
println(df_sol)


df_algo, df_sol, Sol = main(true,wrt=:λ, Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)))
#println(df_algo)
println(df_sol)
df_algo, df_sol, Sol = main(false,wrt=:λ, Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)))
#println(df_algo)
println(df_sol)


df_algo, df_sol, Sol = main(true,wrt=:t0)#,Algorithmes = ((RK4(),3),))
#df_algo, df_sol, Sol = main(true,wrt=:t0,Algorithmes = ((RK4(), 3),))# (Tsit5(), 5),(RadauIIA5(), 9)))
#println(df_algo)
println(df_sol)
df_algo, df_sol, Sol = main(false,wrt=:t0)#,Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)))
#println(df_algo)
println(df_sol)

df_algo, df_sol, Sol = main(true,wrt=:tf)#, Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)))
#println(df_algo)
println(df_sol)
df_algo, df_sol, Sol = main(false,wrt=:tf)
#println(df_algo)
println(df_sol)
=#