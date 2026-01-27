# Examples of second member
#
# Brusselator
# ref : Hairer tome 1 page 201
#=
tf = 20
tspan = (0.0, tf)
x₁0 = 1.3
par = [3.0]
tol = 1.e-4
=#
function bruss(x, par, t)
    λ = par[1]
    x₁ = x[1]
    x₂ = x[2]
    return [1+x₁^2*x₂-(λ+1)*x₁, λ*x₁-x₁^2*x₂]
end

# Linear system
# -------------
# Derivative of the flow example 1
# ẋ = [λ₁x₁ ; λ₂x₂ ; (λ₁-λ₂)x₃]
# x(t0) = x₀(λ) = [λ₂, 1, 1]
# x(t,t0,x_0(λ),λ) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))x₀(λ)
# (∂x(t,t0,x_0(λ),λ)/∂x0) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))
# ∂x(t,t0,x_0(λ),λ)/∂λ = [λ₂texp(λ₁t) exp(λ₁t) ; 0 texp(λ₂t) ; texp((λ₁-λ₂)t)) -texp((λ₁-λ₂)t))]
# 
t0 = 0.
A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
fun_lin(x,λ,t) = A(λ)*x
# Initial condition
fun_lin_x0(λ) = [λ[2],1.,1]; 
# Derivatives of the flow
sol_∂xO_flow(tf,λ) = sol_∂xO_flow = exp((tf-t0)*A(λ))

tf = 2.

sol_∂λ_flow(tf,λ) = [λ[2]*exp(λ[1]*tf) exp(λ[1]*tf)
                    0.           tf*exp(λ[2]*tf)
             tf*exp((λ[1]-λ[2])*tf)  -tf*exp((λ[1]-λ[2])*tf)]

#= 
λ = [1.,2]; p = length(λ);
t0 = 0. ; tf = 1.; tspan = (t0,tf);
funx0(λ) = [λ[2],1.,1]; 
n = length(funx0(λ))
=#


