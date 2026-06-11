using Pkg
Pkg.activate(".")
#Pkg.activate("../")
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

using OrdinaryDiffEq


#include("./fun_examples.jl")
include("../src/CTDiffFlow.jl")
using .CTDiffFlow
include("../src/myode43/myode43.jl")

function test_FD!(df_sol, fun::Function, tspan::Tuple{<:Real,<:Real}, x0::Vector{<:Real}, λ::Vector{<:Real},
         sol_∂xO_flow::Matrix{<:Real},
         adaptive::Bool;
         internalnorm = :default,
         Algorithmes = (("myode43", 7), (RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9)),
         VarInd = ("IND", "VAR2", "VAR1"),
         reltol = 1.e-3, abstol = 1.e-6)
  """
    Tests of different algorithmes of numerical integration 
    the backend ForwardDiff for the automatic differentiation
    input
    -----
    fun, tspan, x0 and λ : definition of the ode
    adaptive : Boolean
               true  : variable steps
               false : fixed steps
    internalnorm = function which computes the norm for controlling the step
    Algorithmes : Tuples
                  Algorithmes[i] : Tuples = name of an algorithm, number of caracters for print
    reltol and abstol : classical relative and absolute error
    return
    ------
    df_odedf_sol : dataframe
    Sol : 
  """
  
    tol_error = 2*eps()
    tol_error = 10*max(reltol,abstol)
    t0 = tspan[1]; tf = tspan[2];
    dt = (tf-t0)/20
   # df_sol = DataFrame(adaptive=Bool[], VAR_IND=String[], internalnorm=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[])
    Sol = []
# Test of automatic differentiation
# ---------------------------------
    ind_Sol = 1
    for algorithme in Algorithmes

      algo = algorithme[1]
      name_algo = string(algorithme[1])[1:algorithme[2]]
    # Integration of the IVP
      if algo == "myode43"
        if adaptive
            T,X = myode43(fun,x0,λ,(t0,tf),reltol,abstol)
        else
            T,X = myode43(fun,x0,λ,t0:dt:tf)
        end  
        push!(df_sol, [adaptive, name_algo, string(Symbol(internalnorm)), NaN, NaN, T[2:3]])
      else
        ivp = ODEProblem(fun, x0, (t0,tf), λ)
        if adaptive
            sol = solve(ivp, alg=algo, reltol = reltol, abstol = abstol)
        else
            sol = solve(ivp, alg=algo, adaptive=adaptive, dt=dt)
        end
        push!(df_sol, [adaptive, name_algo, string(Symbol(internalnorm)), NaN, NaN, sol.t[2:3]])
        T = sol.t
      end

      ind = 1
      for var_ind in VarInd
        if var_ind == "IND"
          RelTol = reltol
          AbsTol = abstol
          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ; var_ind=:ind)
        elseif var_ind =="VAR2"
          RelTol = reltol
          AbsTol = abstol
          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ)
        else
          n = length(x0)
          RelTol = reltol*ones(n,n+1)
          AbsTol = abstol*ones(n,n+1)
          ∂x0_flow = CTDiffFlow.build_∂x0_flow(fun, t0, x0, tf, λ)
        end
        if adaptive
            #sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, reltol=RelTol, abstol=AbsTol)
            

            if internalnorm == :default
              sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, adaptive=true, reltol=RelTol, abstol=AbsTol)
            else
              sol, T = ∂x0_flow(t0, x0, tf, λ; internalnorm = internalnorm, print_times=true, alg=algo, adaptive=true, reltol=RelTol, abstol=AbsTol)
            end
        else
            sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, adaptive=adaptive, dt=dt, alg=algo, reltol=RelTol, abstol=AbsTol)
        end
        push!(Sol,sol)

        norm_inf = norm(sol-sol_∂xO_flow,Inf)
        if ind==1
            norm_diff = NaN
            ind_Sol = ind_Sol+1
        else
            norm_diff = norm(Sol[ind_Sol]-Sol[ind_Sol-1],Inf)
            ind_Sol = ind_Sol+1
        end
        ind = ind+1
        push!(df_sol, [adaptive, var_ind, string(Symbol(internalnorm)), norm_inf, norm_diff, T[2:3]])
      #println(@test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))
       #  @test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error)
      end
    end
    return df_sol,Sol
  end

#using SciMLSensitivity

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


df_sol = DataFrame(adaptive=Bool[], VAR_IND=String[], internalnorm=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[])

test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true)

println(df_sol)

# with my_norm the diagram switches 
sse(x::Number) = x^2
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) #+ sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
function totallength(x::ForwardDiff.Dual)
  totallength(ForwardDiff.value(x)) #+ sum(totallength, ForwardDiff.partials(x))
end
totallength(x::AbstractArray) = sum(totallength, x)
function my_norm(u, t)
  return sqrt(sum(x -> sse(x), u) / totallength(u))
end

#df_sol, Sol = test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=my_norm)#, Algorithmes = ((RK4(), 3),))
#println(df_sol)

#=

# in the automatic differentiation of the flow there is h'(p) the step derivative, 
# so the diagram doesn't switch

df_sol = DataFrame(adaptive=Bool[], VAR_IND=String[], internalnorm=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[])
#algo = (("myode43", 7),(RK4(), 3), (Tsit5(), 5), (RadauIIA5(), 9))
algo = ((RadauIIA5(), 9),)

#df_sol, Sol = test_FD(fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=:default)#,Algorithmes = ((RK4(), 3),))
#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=(u,t) -> norm(u)/sqrt(length(u)),Algorithmes = algo)
#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=(u,t) -> norm(u)/sqrt(length(u)),Algorithmes = ((RK4(), 3),))

#
test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=:default,Algorithmes = algo, VarInd = ("VAR2", "VAR1"))

var_ind = ("IND", "VAR2")
test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=:default,Algorithmes = algo, VarInd = var_ind)
#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=:default,Algorithmes = ((RK4(), 3),))

#println(df_sol)
#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,false,internalnorm=:default)#,Algorithmes = ((RK4(), 3),))
#println(df_sol)

# my_norm2 is the norm of 
# https://github.com/ODINN-SciML/DiffEqSensitivity-Review/blob/main/code/SensitivityForwardAD/example-AD-tolerances.jl

sse2(x::Number) = x^2
sse2(x::ForwardDiff.Dual) = sse2(ForwardDiff.value(x)) + sum(sse2, ForwardDiff.partials(x))
totallength2(x::Number) = 1
function totallength2(x::ForwardDiff.Dual)
  totallength2(ForwardDiff.value(x)) + sum(totallength2, ForwardDiff.partials(x))
end
totallength2(x::AbstractArray) = sum(totallength2, x)
function my_norm2(u, t)
  return sqrt(sum(x -> sse2(x), u) / totallength2(u))
end
test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=my_norm2, Algorithmes = algo, VarInd = var_ind)
#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=my_norm2, Algorithmes = ((RK4(), 3),))
#


#

test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=my_norm, Algorithmes = algo, VarInd = var_ind)

#test_FD!(df_sol,fun_lin, tspan, x0, λ, sol_∂xO_flow,true,internalnorm=my_norm, Algorithmes = ((RK4(), 3),))


println(df_sol)




=#