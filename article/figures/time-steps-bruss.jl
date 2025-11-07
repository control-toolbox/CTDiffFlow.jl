
using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings
using ForwardDiff: ForwardDiff
using Zygote: Zygote
#using Mooncake # plante à l'instalation
using SciMLSensitivity
using DifferentiationInterface

include("../../src/CTDiffFlow.jl")
using .CTDiffFlow
#
# Definition of the second member
function bruss(x,par,t)
  λ = par[1]
  x₁ = x[1]
  x₂ = x[2]
  return [1+x₁^2*x₂-(λ+1)*x₁ , λ*x₁-x₁^2*x₂]
end

t0 = 0.; tf = 1.
tspan = (t0,tf)
λ = [3.]
p = length(λ)
x01 = 1.3
funx0(λ) = [x01, λ[1]]
tol = 1.e-4
n = 2

plt1 = plot(); plt2 = plot()
algo = Tsit5()
reltol = 1.e-4; abstol = 1.e-4

# Flow
ivp = ODEProblem(bruss, funx0(λ), (t0,tf), λ)
sol_ivp = solve(ivp, alg=algo; reltol=reltol,abstol=abstol)
println("Flow")
println("Times for the initial flow T = ")
println(sol_ivp.t)
dict_T = Dict(:ivp => sol_ivp.t)
dt = sol_ivp.t[2]

# -------------------------------
println("Vatiational equation")
# ------------------------------

∂λ_flow_var = CTDiffFlow.build_∂λ_flow_var(bruss,t0,funx0,tf,λ)

#println(∂λ_flow_var(t0,funx0,tf,λ;reltol=reltol,abstol=abstol))
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=reltol,abstol=abstol)
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=reltol,abstol=abstol,internalnorm = (u,t)->norm(u))
println("Times default values= ", sol_var.t)
dict_T[:var_reltol] = sol_var.t
my_Inf = prevfloat(typemax(Float64))
RelTol = [reltol*ones(n,1) Inf*ones(n,1)]/sqrt(p+1)
AbsTol = [abstol*ones(n,1) Inf*ones(n,1)]/sqrt(p+1)
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=RelTol,abstol=AbsTol)
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=RelTol,abstol=AbsTol, dt=dt)#,internalnorm = (u,t)->norm(u) )
println("RelTol = ", RelTol)
println("Times = ", sol_var.t)
println(sol_ivp.t - sol_var.t)
dict_T[:var_RelTol] = sol_var.t

# ne fonctionne pas 
#@btime ∂λ_flow_var(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

println("Automatic differentiation")
println("with ForwardDiff")
∂λ_flow = CTDiffFlow.build_∂λ_flow(bruss,t0,funx0,tf,λ)
sol_diff_auto_flow, T = ∂λ_flow(t0,funx0,tf,λ;reltol=RelTol,abstol=AbsTol,print_times=true)
println("Times : ", T)
dict_T[:diff_auto_ForwardDiff] = T
#sol_diff_auto_flow, T = ∂λ_flow(tspan,funx0,λ;reltol=RelTol,abstol=AbsTol,print_times=true)
#= 
println("internalnorm = (u,t)->norm(u)")
println(∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,internalnorm = (u,t)->norm(u),print_times=true))
#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)
=#
println("With Zygote")
# with Zygote
∂λ_flow = CTDiffFlow.build_∂λ_flow(bruss,t0,funx0,tf,λ; backend=AutoZygote())
sol_diff_auto_flow_tf, T = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,print_times=true)
dict_T[:diff_auto_Zygote] = T
println(sol_diff_auto_flow_tf-sol_var)
#sol_diff_auto_flow,T = ∂λ_flow(tspan,funx0,λ;reltol=reltol,abstol=abstol,print_times=true)
#println(sol_diff_auto_flow)
#reshape(sol_diff_auto_flow,2,7)

#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

println("With Mooncake")
# plante
#∂λ_flow = CTDiffFlow.build_∂λ_flow(bruss,t0,funx0,tf,λ; backend=AutoMooncake())
#println(∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,print_times=true))
#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

#plot(plt1,plt2)

