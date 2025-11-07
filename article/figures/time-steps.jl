
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

# Derivative of the flow example 1
# ẋ = [λ₁x₁ ; λ₂x₂ ; (λ₁-λ₂)x₃
# x(t0) = x₀(λ) = [λ₂, 1, 1]
# x(t,t0,x_0(λ),λ) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))x₀(λ)
# (∂x(t,t0,x_0(λ),λ)/∂x0) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))
# ∂x(t,t0,x_0(λ),λ)/∂λ = [λ₂texp(λ₁t) exp(λ₁t) ; 0 texp(λ₂t) ; texp((λ₁-λ₂)t)) -texp((λ₁-λ₂)t))]
# 
# ∂x0_flow

A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
fun1(x,λ,t) = A(λ)*x
function my_∂x0_flow(tspan,x0,λ; ode_kwargs...)
    n = length(x0)
    In = Matrix(I(n))
    ivp_fun1 = ODEProblem(fun1, In , tspan, λ)
    alg = get(ode_kwargs, :alg, Tsit5())
    sol = solve(ivp_fun1, alg=alg; ode_kwargs...)
    return sol
end

function my_∂x0_flow(t0,x0,tf,λ; ode_kwargs...)
    sol = my_∂x0_flow((t0,tf),x0,λ; ode_kwargs...)
    return sol.u[end]
end

function my_∂λ_flow(tspan,xX0,λ; ode_kwargs...)
    #n,p = size(xX0)
    fun2(xX,λ,t) = [fun1(xX[:,1],λ,t) fun1(xX[:,2:end],λ,t)+[xX[1,1] 0 ; 0 xX[2,1] ; xX[3,1] xX[3,1]] ]
    
    ivp_fun = ODEProblem(fun1, xX0 , tspan, λ)
    alg = get(ode_kwargs, :alg, Tsit5())
    sol = solve(ivp_fun, alg=alg; ode_kwargs...)
    return sol
end

function my_∂λ_flow(t0,xX0,tf,λ; ode_kwargs...)
    sol = my_∂λ_flow((t0,tf),xX0,λ; ode_kwargs...)
    return sol.u[end]
end


λ = [1.,2]; p = length(λ);
t0 = 0. ; tf = 1.; tspan = (t0,tf);
funx0(λ) = [λ[2],1.,1]; 
n = length(funx0(λ))
sol_∂xO_flow = diagm([exp(tf*λ[1]),exp(tf*λ[2]),exp(tf*(λ[1]-λ[2]))])
sol_∂λ_flow = [λ[2]*exp(λ[1]*tf) exp(λ[1]*tf)
                    0.           tf*exp(λ[2]*tf)
             tf*exp((λ[1]-λ[2])*tf)  -tf*exp((λ[1]-λ[2])*tf)]
algo = Tsit5()
reltol = 1.e-4; abstol = 1.e-4

# Flow
#x0λ = funx0(λ)
ivp = ODEProblem(fun1, funx0(λ), (t0,tf), λ)
sol_ivp = solve(ivp, alg=algo; reltol=reltol,abstol=abstol)
println("Flow")
println("Times for the initial flow T = ")
println(sol_ivp.t)
dict_T_var = Dict(:ivp => sol_ivp.t)
dict_T_auto_diff = Dict(:ivp => sol_ivp.t)
dt = sol_ivp.t[2]

# -------------------------------
println("Vatiational equation")
# ------------------------------

xX0 = [funx0(λ) [0 1 ; 0 0 ; 0 0]]
println(my_∂λ_flow(t0,xX0,tf,λ))

RelTol = [reltol*ones(n,1) Inf*ones(n,2)]/sqrt(p+1)
AbsTol = [abstol*ones(n,1) Inf*ones(n,2)]/sqrt(p+1)

# plante car 0.*Inf = NaN lors de l'initialisation du pas
sol1 = my_∂λ_flow((t0,tf),xX0,λ;reltol=RelTol,abstol=AbsTol, dt=dt)

∂λ_flow_var = CTDiffFlow.build_∂λ_flow_var(fun1,t0,funx0,tf,λ)

sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=reltol,abstol=abstol)
dict_T_var[:var_reltol] = sol_var.t
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=reltol,abstol=abstol,internalnorm = (u,t)->norm(u[:,1:1]))
dict_T_var[:var_reltol_internalnorm] = sol_var.t

my_Inf = prevfloat(typemax(Float64))
RelTol = [reltol*ones(n,1) my_Inf*ones(n,2)]/sqrt(p+1)
AbsTol = [abstol*ones(n,1) my_Inf*ones(n,2)]/sqrt(p+1)


println("AbsTol = ", AbsTol)
println("RelTol = ", RelTol)
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=RelTol,abstol=AbsTol)
println(sol_ivp.t - sol_var.t)
#println("dt = ",dt)
sol_var = ∂λ_flow_var((t0,tf),funx0,λ;reltol=RelTol,abstol=AbsTol,dt=dt)
println("RelTol = ", RelTol)
println("Times = ", sol_var.t)
println(sol_ivp.t - sol_var.t)
dict_T_var[:var_RelTol] = sol_var.t
xX_var_tf = sol_var.u[end][:,2:end]


# ne fonctionne pas 
#@btime ∂λ_flow_var(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

println("Automatic differentiation")
println("with ForwardDiff")
∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1,t0,funx0,tf,λ)
sol_diff_auto_flow, T = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,print_times=true)
println("Times : ", T)
dict_T_auto_diff[:diff_auto_ForwardDiff] = T

#sol_diff_auto_flow, T = ∂λ_flow(tspan,funx0,λ;reltol=reltol,abstol=abstol,print_times=true)

println("internalnorm = (u,t)->norm(u)")
sol_diff_auto_flow, T = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,internalnorm = (u,t)->norm(u),print_times=true)
println("Times : ", T)
dict_T_auto_diff[:diff_auto_ForwardDiff_internalnorm] = T
#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

println("With Zygote")
# with Zygote
∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1,t0,funx0,tf,λ; backend=AutoZygote())
println("t0= ", t0)
println("tf= ", tf)
println("λ = ", λ)
N = 5
sol_diff_auto_flow_tf, T = ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol, print_times=true) #, alg = Euler(),adaptive=false, dt = (tf-t0)/N,print_times=true)
dict_T_auto_diff[:diff_auto_Zygote] = T

[sol_∂λ_flow sol_diff_auto_flow sol_diff_auto_flow_tf]
#sort(dict_T; rev = true)

#sol_diff_auto_flow,T = ∂λ_flow(tspan,funx0,λ;reltol=reltol,abstol=abstol,print_times=true)
#println(sol_diff_auto_flow)
#reshape(sol_diff_auto_flow,2,7)

#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)

println("With Mooncake")
# plante
#∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1,t0,funx0,tf,λ; backend=AutoMooncake())
#println(∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol,print_times=true))
#@btime ∂λ_flow(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)



