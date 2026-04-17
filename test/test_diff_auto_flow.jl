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

function main(adaptive,internalnorm = :default)
    tol_error = 2*eps()
    reltol = 1.e-3;
    abstol = 1.e-6
    tol_error = 10*max(reltol,abstol)
    #
    # Initial value problem
    λ = [1.0, 2]
    B(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
    fun_lin2(x,λ,t) = B(λ)*x
    t0 = 0. ; tf = 1.
    x0 = [1., 2., 3]
    # jacobien of the flow
    sol_∂xO_flow2(tf,λ) = exp(tf*B(λ))


    # Test of convergence with different the numerical integration algorithms
    Algorithmes = ((RK4(), 3), (Tsit5(), 5),(RadauIIA5(), 9))

    df_algo = DataFrame(VAR_IND=String[], backend=String[], adaptive=Bool[], cv=Bool[])
    

    Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
    #Backends = (AutoForwardDiff(), AutoMooncake(), AutoZygote())
    for algorithme in Algorithmes
        algo = algorithme[1]
        name_algo = string(algorithme[1])[1:algorithme[2]]
        # Integration of the IVP
        if algo == "myode43"
        T,X = myode43(fun_lin2,x0,λ,(t0,tf),reltol,abstol)
        push!(df_algo, ["IVP", name_algo, true])
        else
        ivp = ODEProblem(fun_lin2, x0, (t0,tf), λ)
        sol = solve(ivp, alg=algo, reltol = reltol, abstol = abstol)
        push!(df_algo, ["IVP", name_algo, adaptive, true])
        T = sol.t
        end
        for var_ind in ("IND", "VAR2", "VAR1")
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
              try
                  sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, adaptive=adaptive, reltol=RelTol, abstol=AbsTol)
                  push!(df_algo, [var_ind, string(backend), adaptive, true])
              catch
                  push!(df_algo, [var_ind, string(backend), adaptive, false])
              end
          end
        end
    end
      

# Test of automatic differentiation
# ---------------------------------
    
    # Backends
    # Problems with Enzyme and Zygote
    #Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
    Backends = (AutoForwardDiff(), AutoZygote())
  
    #Backends = (AutoMooncake(),AutoZygote())
    # algorithme = (RK4(), 3)
    algorithme = (Tsit5(), 5)
    #algorithme = (RadauIIA5(), 9)
    # algorithme = ("myode43", 7)
    algo = algorithme[1]
    name_algo = string(algorithme[1])[1:algorithme[2]]
    
    df_sol = DataFrame(VAR_IND=String[], backend=String[], norm_∞_error=Real[], norm_∞_diff=Real[], time_steps=Vector[])
    Sol = []
    # Integration of the IVP
    if algo == "myode43"
      T,X = myode43(fun_lin2,x0,λ,(t0,tf),reltol,abstol)
      println("T = ", T)
      push!(df_sol, ["IVP", name_algo, NaN, NaN, T[2:3]])
    else
      ivp = ODEProblem(fun_lin2, x0, (t0,tf), λ)
      sol = solve(ivp, alg=algo, reltol = reltol, abstol = abstol)
      push!(df_sol, ["IVP", name_algo, NaN, NaN, sol.t[2:3]])
      T = sol.t
    end
    dt = T


      ind = 1
      for var_ind in ("IND", "VAR2", "VAR1")
      #for var_ind in ("IND",)
        for backend in Backends
          println("bachend = ", backend)
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
            #sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, reltol=RelTol, abstol=AbsTol)
            

            if internalnorm == :default
              sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, alg=algo, adaptive=true, reltol=RelTol, abstol=AbsTol)
            else
              sol, T = ∂x0_flow(t0, x0, tf, λ; internalnorm = internalnorm, print_times=true, alg=algo, adaptive=true, reltol=RelTol, abstol=AbsTol)
            end
        else
          println("adaptive = ", adaptive)
            sol, T = ∂x0_flow(t0, x0, tf, λ; print_times=true, adaptive=adaptive, alg=algo, reltol=RelTol, abstol=AbsTol)
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
        push!(df_sol, [var_ind, string(backend), norm_inf, norm_diff, T[2:3]])
      #println("∂x0_flow = ", ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
      #println("sol_∂xO_flow(tf,λ) = ", sol_∂xO_flow(tf,λ))
      #println(@test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))
       #  @test isapprox(sol_∂xO_flow(tf,λ), ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error)
      end
    end

    return df_algo, df_sol,Sol
  end

using SciMLSensitivity

# in the automatic differentiation of the flow there is h'(p) the step derivative, 
# so the diagram doesn't switch

df_algo, df_sol, Sol = main(true)
#println(df_algo)
println(df_sol)

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
df_algo, df_sol, Sol = main(true, my_norm2)
#println(df_algo)
println(df_sol)

# with my_norm the diagram switches 
sse(x::Number) = x^2
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) #+ sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
function totallength(x::ForwardDiff.Dual)
  totallength(ForwardDiff.value(x)) #+ sum(totallength, ForwardDiff.partials(x))
end
totallength(x::AbstractArray) = sum(totallength, x)
my_norm = (u, t) -> sqrt(sum(x -> sse(x), u) / totallength(u))

df_algo, df_sol, Sol = main(true, my_norm)
#println(df_algo)
println(df_sol)




# the diagram switches
df_algo, df_sol, Sol = main(false)
println(df_algo)
println(df_sol)



