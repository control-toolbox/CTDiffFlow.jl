using Pkg
Pkg.activate(".")
#Pkg.add("BenchmarkTools")
#Pkg.add("Enzyme")
#Pkg.add("Mooncake")
#println(pwd())
#println(Pkg.status())
using Markdown
using LinearAlgebra
using Test
using DifferentiationInterface
using ForwardDiff: ForwardDiff
using Enzyme: Enzyme
using Mooncake: Mooncake
using Zygote: Zygote

function main()
A = [0. 1 2 ; 3 4 5]
function fun1(x::Vector{<:Real})::Vector{<:Real}
    return A*x
end
Jac_fun1 = A

B = [0. 1 2 ; 3 4 5 ; 6 7 8]
function fun2(x::Vector{<:Real})::Vector{<:Real}
    return [0.5*transpose(x)*B*x]
end
B_sym = 0.5*(transpose(B)+B)

function fun3(x::Matrix{<:Real})::Matrix{<:Real}
    return A*x
end
Jac_fun3 = [A zeros(size(A)) ; zeros(size(A)) A]

tol_error = 2*eps()

# Vectors
println("Vectors")
println("-------")
# Problems with Enzyme
#Backends = (AutoEnzyme(), AutoForwardDiff(), AutoMooncake(), AutoZygote())
Backends = (AutoForwardDiff(), AutoMooncake(), AutoZygote())
    println("type_x = Vector")
    for backend in Backends
       println("backend = ", backend)
       x = [1.,2,3]
       Jac = jacobian(fun1, backend, x)
       println(@test isapprox(Jac_fun1, Jac, atol=tol_error))
       Jac = jacobian(fun2, backend, x)
       println(@test isapprox(reshape(B_sym*x,1,3), Jac, atol=tol_error))
    end

# Matrix
println("Matrix")
println("------")
# Problem with Mooncake
Backends = (AutoForwardDiff(), AutoZygote())
for backend in Backends
    println("backend = ", backend)
    x = [1 2 ; 3 4 ; 5 6]
    Jac = jacobian(fun3, backend, x)
    #println("true_jac = ", Jac_fun3)
    #println("Jac_autodiff = ", Jac)
    println(@test isapprox(Jac_fun3, Jac, atol=tol_error))
end


#=

# benchmark
using BenchmarkTools

println(@benchmark jacobian($fun1, $backend1, $x))
# with preparation
prep1 = prepare_jacobian(fun1, backend1, x)
Jac_jun1 = similar(A)
jacobian!(fun1, Jac_fun1, prep1, backend1, x)
Jac_fun1  # has been mutated

=#

end
main()
