
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using Zygote: Zygote
#include("DiffFlow.jl")
#using .DiffFlow
#using OrdinaryDiffEq

function tt(a::Vector{<:Real}; b=1.0::Real)
    return a[1]+b
end

B = [1 2; 3 4]
fun1(x) = B*x
t0 = 0.0;
tf = 1.0;
x0 = [1, 2];
tspan = (t0, tf)
backend=AutoZygote()
backend1 = AutoForwardDiff()
jacobian(fun1, backend1, x0)
