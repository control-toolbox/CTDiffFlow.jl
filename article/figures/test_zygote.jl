using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings
using ForwardDiff: ForwardDiff
using Zygote: Zygote
#using Mooncake # plante à l'instalation
using SciMLSensitivity
using DifferentiationInterface
using DataFrames

include("../../src/CTDiffFlow.jl")
using .CTDiffFlow
include("../../src/myode43/myode43.jl")

function my_euler(fun,t0,x0,tf,λ,N)
    ti = t0; xi = x0
    h = (tf-t0)/N
    for i in 1:N
        xi +=  h*fun(xi,λ,ti)
        ti = ti+h
    end
    return xi
end


# Definition of the edo
A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
fun1(x,λ,t) = A(λ)*x
funx0(λ) = [λ[2],1.,1]; 
N = 10
t0 = 0. ; tf = 1.; tspan = (t0,tf);
λ = [1.,2];
reltol = 1.e-8
abstol = reltol
# ∂λ_flow(λ)
sol_∂λ_flow = [λ[2]*exp(λ[1]*tf) exp(λ[1]*tf)
                    0.           tf*exp(λ[2]*tf)
             tf*exp((λ[1]-λ[2])*tf)  -tf*exp((λ[1]-λ[2])*tf)]


#dict_Zygote = Dict(:∂λ_sol => sol_∂λ_flow)
#dict_ForwardDiff = Dict(:∂λ_sol => sol_∂λ_flow)

df = DataFrame(AutoDiff = String[], algo=String[], Jacobian = Matrix{<:Real}[])
push!(df, ("solution", "", sol_∂λ_flow))

# Zygote with my_euler
_flow(λ) = my_euler(fun1,t0,funx0(λ),tf,λ,N)
#dict_Zygote[:my_euler] = jacobian(_flow,AutoZygote(),λ)
#dict_ForwardDiff[:my_euler] = jacobian(_flow,AutoForwardDiff(),λ)
push!(df, ("Zygote", "my-euler", jacobian(_flow,AutoZygote(),λ)))
push!(df, ("ForwardDiff", "my-euler", jacobian(_flow,AutoForwardDiff(),λ)))
# Zygote with numerical integration
function _flow_int(λ)
            ivp = ODEProblem(fun1, funx0(λ), (t0,tf), λ)
            #algo = get(ode_kwargs, :alg, Tsit5())
            #println("algo = ", algo)

            sol = solve(ivp, alg = Euler(),reltol=reltol,abstol=abstol,adaptive=false, dt = (tf-t0)/N)
            return sol.u[end]
end

#dict_Zygote[:Euler] = jacobian(_flow_int,AutoZygote(),λ)
#dict_ForwardDiff[:Euler] = jacobian(_flow_int,AutoForwardDiff(),λ)
push!(df, ("Zygote", "Euler", jacobian(_flow_int,AutoZygote(),λ)))
push!(df, ("ForwardDiff", "Euler", jacobian(_flow_int,AutoForwardDiff(),λ)))
# Zygote with ∂λ_flow
∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1,t0,funx0,tf,λ; backend=AutoZygote())
#dict_Zygote[:CTDiffFlowZygoteEuler] = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol, alg = Euler(),adaptive=false, dt = (tf-t0)/N)#,print_times=false)
push!(df, ("Zygote", "CTDiffFlow-Euler", ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol, alg = Euler(),adaptive=false, dt = (tf-t0)/N)))

# ForwardDiff with ∂λ_flow
∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1,t0,funx0,tf,λ; backend=AutoForwardDiff())
#dict_ForwardDiff[:CTDiffFlowForwardDiffEuler] = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol, alg = Euler(),adaptive=false, dt = (tf-t0)/N)#,print_times=false)
push!(df, ("ForwardDiff", "my-euler", ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol, alg = Euler(),adaptive=false, dt = (tf-t0)/N)))

function _flow(λ)
    T, X = myode43(fun1,funx0(λ),λ,(t0,tf),abstol,reltol)
    return X[end]
end

#dict_Zygote[:myode43Zygote] = jacobian(_flow,AutoZygote(),λ)
#dict_ForwardDiff[:myode43ForwardDiff] = jacobian(_flow,AutoForwardDiff(),λ)
push!(df, ("Zygote", "my-ode43", [NaN;;]))
push!(df, ("ForwardDiff", "my-ode43", jacobian(_flow,AutoForwardDiff(),λ)))

jacobian(_flow,AutoForwardDiff(),λ)

latexify(df,arraystyle=:pmatrix,env=:tabular)

sol_Zygote = df[in.(df.AutoDiff, Ref(["solution","Zygote"])),:]
sol_ForwardDiff = df[in.(df.AutoDiff, Ref(["solution", "ForwardDiff"])),:]
