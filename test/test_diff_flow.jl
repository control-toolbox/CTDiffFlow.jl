using Pkg
Pkg.activate(".")
#Pkg.add("DualNumbers")
using DualNumbers
using Markdown
using LinearAlgebra
using Test
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using Zygote: Zygote
using OrdinaryDiffEq

include("../src/CTDiffFlow.jl")
using .CTDiffFlow

function rhs_bruss(x, par, t)
    λ = par[1]
    x₁ = x[1]
    x₂ = x[2]
    return [1+x₁^2*x₂-(λ+1)*x₁, λ*x₁-x₁^2*x₂]
end

function my_rhs_bruss_var_x(xX, Λ, t)
    λ = Λ[1]
    x = xX[:, 1]
    x₁ = x[1]
    x₂ = x[2]
    X = xX[:, 2:end]
    xpoint = rhs_bruss(x, Λ, t)

    rhs_bruss_x = [
        2*x₁*x₂-λ-1 x₁^2
        λ-2*x₁*x₂ -x₁^2
    ]
    Xpoint = rhs_bruss_x*X

    return [xpoint Xpoint]
end

function my_rhs_bruss_var_λ(xX, Λ, t)
    λ = Λ[1]
    x = xX[:, 1]
    x₁ = x[1]
    x₂ = x[2]
    X = xX[:, 2:end]
    xpoint = rhs_bruss(x, Λ, t)

    rhs_bruss_x = [
        2*x₁*x₂-λ x₁^2
        λ-2*x₁*x₂ -x₁^2
    ]
    rhs_bruss_λ = [-x₁; x₁]
    Xpoint = rhs_bruss_x*X + rhs_bruss_λ

    return [xpoint Xpoint]
end

tol_error = 1.e-5

# Tests of the second member of the variational equations
# with respect to the initial condition
rhs_bruss_var_x = CTDiffFlow.built_rhs_var(rhs_bruss)
t = 0.0
Λ = [3.0]
x0 = [1.3; Λ[1]]
xX = [x0 [0; 1]]
@test isapprox(my_rhs_bruss_var_x(xX, Λ, t), rhs_bruss_var_x(xX, Λ, t), atol=tol_error)
xX = [x0 [1; 0]]
@test isapprox(my_rhs_bruss_var_x(xX, Λ, t), rhs_bruss_var_x(xX, Λ, t), atol=tol_error)

# with respect to λ
rhs_bruss_var_λ = CTDiffFlow.built_rhs_var(rhs_bruss; wrt=:λ)
t = 0.0
Λ = [3.0]
x0 = [1.3; Λ[1]]
xX = [x0 [0; 1]]
@test isapprox(my_rhs_bruss_var_λ(xX, Λ, t), rhs_bruss_var_λ(xX, Λ, t), atol=tol_error)
#xX = [x0 [1;0]]
#@test isapprox(my_rhs_bruss_var_λ(xX,Λ,t), rhs_bruss_var_λ(xX,Λ,t) , atol = tol_error) 

# Derivative of the flow example 1
# ẋ = [λ₁x₁ ; λ₂x₂ ; (λ₁-λ₂)x₃
# x(t0) = x₀(λ) = [λ₂, 1, 1]
# x(t,t0,x_0(λ),λ) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))x₀(λ)
# (∂x(t,t0,x_0(λ),λ)/∂x0) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))
# ∂x(t,t0,x_0(λ),λ)/∂λ = [λ₂texp(λ₁t) exp(λ₁t) ; 0 texp(λ₂t) ; texp((λ₁-λ₂)t)) -texp((λ₁-λ₂)t))]
# 
# ∂x0_flow
A(λ) = [λ[1] 0 0; 0 λ[2] 0; 0 0 λ[1]-λ[2]]
fun1(x, λ, t) = A(λ)*x
reltol = 1.e-8;
abstol = 1.e-12
function my_∂x0_flow(tspan, x0, λ; ode_kwargs...)
    n = length(x0)
    In = Matrix(I(n))
    ivp_fun1 = ODEProblem(fun1, In, tspan, λ)
    alg = get(ode_kwargs, :alg, Tsit5())
    sol = solve(ivp_fun1; alg=alg, ode_kwargs...)
    return sol
end

function my_∂x0_flow(t0, x0, tf, λ; ode_kwargs...)
    sol = my_∂x0_flow((t0, tf), x0, λ; ode_kwargs...)
    return sol.u[end]
end

λ = [1.0, 2]
t0 = 0.0;
tf = 2.0;
x0 = [λ[2], 1.0, 1];
sol_∂xO_flow = exp((tf-t0)*A(λ))
# sol_∂xO_flow = diagm([exp(tf*λ[1]), exp(tf*λ[2]), exp(tf*(λ[1]-λ[2]))])

∂x0_flow_var = CTDiffFlow.build_∂x0_flow_var(fun1, t0, x0, tf, λ)
println(@test isapprox(sol_∂xO_flow, ∂x0_flow_var(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))

println(my_∂x0_flow(t0, x0, tf, λ))
println(my_∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))

println(my_∂x0_flow(t0, x0, tf, λ) - ∂x0_flow_var(t0, x0, tf, λ))
println(@test isapprox(
    my_∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol),
    ∂x0_flow_var(t0, x0, tf, λ; reltol=reltol, abstol=abstol),
    atol=tol_error))

# Derivative with respect to λ
sol_∂λ_flow = [
    λ[2]*tf*exp(λ[1]*tf) exp(λ[1]*tf)
    0 tf*exp(λ[2]*tf)
    tf*exp((λ[1]-λ[2])*tf) -tf*exp((λ[1]-λ[2])*tf)
]

funx0(λ) = [λ[2], 1.0, 1]
∂λ_flow_var = CTDiffFlow.build_∂λ_flow_var(fun1, t0, funx0, tf, λ)
λ0 = zeros(3, 2);
λ0[1, 2] = 1.0
println(@test isapprox(
    sol_∂λ_flow, ∂λ_flow_var(t0, funx0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error))

# Diff auto
# Derivative with respect to x0
∂x0_flow = CTDiffFlow.build_∂x0_flow(fun1, t0, x0, tf, λ)
println("ccc", sol_∂xO_flow-∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
println(
    "ddd",
    (
        my_∂x0_flow(
            t0, x0, tf, λ; reltol=reltol, abstol=abstol
        )-∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol)
    ),
)
@test isapprox(
    sol_∂xO_flow, ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error
)
@test isapprox(
    my_∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol),
    ∂x0_flow(t0, x0, tf, λ; reltol=reltol, abstol=abstol),
    atol=tol_error,
)

# Derivative with respect to λ

∂λ_flow = CTDiffFlow.build_∂λ_flow(fun1, t0, funx0, tf, λ)
println("comparaison diff auto et var")
println("----------------------------")
println(funx0(λ))
println(∂λ_flow(t0, funx0, tf, λ; print_times=true, reltol=reltol, abstol=abstol))
println(sol_∂λ_flow)
@test isapprox(
    sol_∂λ_flow, ∂λ_flow(t0, funx0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error
)

println(∂λ_flow(t0, funx0, tf, λ; reltol=reltol, abstol=abstol))
println(sol_∂λ_flow)
@test isapprox(
    sol_∂λ_flow, ∂λ_flow(t0, funx0, tf, λ; reltol=reltol, abstol=abstol), atol=tol_error
)
