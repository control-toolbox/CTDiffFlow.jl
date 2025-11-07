include("../../src/myode43/myode43.jl")


# ẋ = [λ₁x₁ ; λ₂x₂ ; (λ₁-λ₂)x₃
# x(t0) = x₀(λ) = [λ₂, 1, 1]
# x(t,t0,x_0(λ),λ) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))x₀(λ)
# (∂x(t,t0,x_0(λ),λ)/∂x0) = diag(exp(λ₁t),exp(λ₂t),exp((λ₁-λ₂)t))
# ∂x(t,t0,x_0(λ),λ)/∂λ = [λ₂texp(λ₁t) exp(λ₁t) ; 0 texp(λ₂t) ; texp((λ₁-λ₂)t)) -texp((λ₁-λ₂)t))]
# 
# ∂x0_flow

A(λ) = [λ[1] 0 0 ; 0 λ[2] 0 ; 0 0 λ[1]-λ[2]]
fun1(x,λ,t) = A(λ)*x

fun2(xX,λ,t) = [fun1(xX[:,1],λ,t) fun1(xX[:,2:end],λ,t)+[xX[1,1] 0 ; 0 xX[2,1] ; xX[3,1] xX[3,1]] ]
   


t0 = 0. ; tf = 2.; tspan = (t0,tf);
λ = [1.,2]; p = length(λ);
funx0(λ) = [λ[2],1.,1];  
x0 = funx0(λ)
n = length(x0)
reltol = 1.e-4; abstol = 1.e-4
myode43(fun1,x0,λ,tspan,reltol,abstol)

xX0 = [x0 [0 1 ; 0 0 ; 0 0]]

myode43(fun2,xX0,λ,tspan,reltol,abstol)

RelTol = reltol*ones(n,p+1) 
AbsTol = RelTol

myode43(fun2,xX0,λ,tspan,RelTol,AbsTol)

RelTol = [reltol*ones(n,1) Inf*ones(n,p)]/sqrt(p+1)
AbsTol = [abstol*ones(n,1) Inf*ones(n,p)]/sqrt(p+1)

inith(fun1,x0,λ,t0,abstol,reltol)

inith(fun2,xX0,λ,t0,abstol,reltol)
RelTol = reltol*ones(n,p+1) 
AbsTol = RelTol
inith(fun2,xX0,λ,t0,AbsTol,RelTol)


epsi = eps()
epsi = 0
xX0 = [x0 [epsi 1 ; epsi epsi ; epsi epsi]]
my_Inf = prevfloat(typemax(Float64))
#my_Inf = Inf
RelTol = [reltol*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
AbsTol = [abstol*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
inith(fun2,xX0,λ,t0,AbsTol,RelTol)
myode43(fun2,xX0,λ,tspan,RelTol,AbsTol)

