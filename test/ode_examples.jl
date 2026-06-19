# Examples of second member
# Linear system
# -------------
# Initial value problem
struct ode_lin
    t0  :: Real
    tf  :: Real
    x0  :: Vector{<:Real}
    λ   :: Vector{<:Real}
    rhs :: Function
    eq_var_x0 :: Function
    eq_var_λ   :: Function
    xf  :: Vector{<:Real}
    # jacobien of the flow at xf
    sol_∂x0_flow :: Matrix{<:Real}
    sol_∂λ_flow  :: Matrix{<:Real}
    sol_∂t0_flow :: Vector{<:Real}
    sol_∂tf_flow :: Vector{<:Real}
end

t0 = 0. ; tf = 1.
x0 = [1., 2., 3]
λ = [1.0, 2]
A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
rhs_lin1(x,λ,t) = A(λ)*x
∂x_rhs_lin1(x,λ,t) = A(λ)
∂λ_rhs_lin1(x,λ,t) = [1 0 0 ; 0 0 0 ; 0 0 1]*x + [0 0 0 ; 0 1 0 ; 0 0 -1]*x
eq_var_x0(xX,λ,t) = [rhs_lin1(xX[:,1],λ,t) A(λ)*xX[:,2:end]]
eq_var_λ(xX,λ,t) = [rhs_lin1(xX[:,1],λ,t) A(λ)*xX[:,2:end]+[[1 0 0 ; 0 0 0 ; 0 0 1]*xX[:,1] [0 0 0 ; 0 1 0 ; 0 0 -1]*xX[:,1]]]
xf = exp((tf-t0)*A(λ))*x0
sol_∂x0_flow = exp((tf-t0)*A(λ))
sol_∂λ_flow = (tf-t0)*[xf[1] 0 ; 0 xf[2] ; xf[3] -xf[3]]
sol_∂t0_flow = [ -xf[1]*λ[1] , -xf[2]*λ[2] , xf[3]*(λ[2]-λ[1])]
sol_∂tf_flow = rhs_lin1(xf,λ,tf)

ode_lin1 = ode_lin(t0,tf,x0,λ,rhs_lin1,eq_var_x0, eq_var_λ, xf,sol_∂x0_flow,sol_∂λ_flow,sol_∂t0_flow,sol_∂tf_flow)


# Brusselator
# ref : Hairer tome 1 page 201



struct bruss_struct
    t0  :: Real
    tf  :: Real
    x0  :: Vector{<:Real}
    λ   :: Vector{<:Real}
    rhs :: Function
    eq_var_x0 :: Function
    eq_var_λ   :: Function
    #=
    xf  :: Vector{<:Real}
    # jacobien of the flow at xf
    sol_x0_∂flow :: Matrix{<:Real}
    sol_λ_∂flow  :: Matrix{<:Real}
    sol_t0_∂flow :: Vector{<:Real}
    sol_tf_∂flow :: Vector{<:Real}
    =#
end
t0 = 0.
tf = 20
λ = [3.0]
x0 = [1.3, λ[1]]
function brusselator(x, par, t)
    λ = par[1]
    x₁ = x[1]
    x₂ = x[2]
    return [1+x₁^2*x₂-(λ+1)*x₁, λ*x₁-x₁^2*x₂]
end

function ∂x_brusselator(x,par,t)
    x₁ = x[1]
    x₂ = x[2]
    λ = par[1]
    return [2*x₁*x₂-(λ+1) x₁^2 ; λ-2*x₁*x₂ -x₁^2]
end
function ∂λ_brusselator(x,par,t)
    x₁ = x[1]
    x₂ = x[2]
    λ = par[1]
    return [-x₁ ; x₁]
end

function eq_var_bruss_x0(xδx,λ,t)
    x = xδx[:,1]
    δx = xδx[:,2:end]
    return [brusselator(x,λ,t) ∂x_brusselator(x,λ,t)*δx]
end
function eq_var_bruss_λ(xδλ,λ,t)
    x = xδλ[:,1]
    δλ = xδλ[:,2:end]
    return [brusselator(x,λ,t) ∂x_brusselator(x,λ,t)*δλ+∂λ_brusselator(x,λ,t)]
end


bruss = bruss_struct(t0,tf,x0,λ,brusselator,eq_var_bruss_x0,eq_var_bruss_λ)