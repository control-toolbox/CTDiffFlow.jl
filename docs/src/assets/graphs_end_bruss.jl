using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK
#using ODEInterfaceDiffEq
using ODEInterface, ODEInterfaceDiffEq
ODEInterface.loadODESolvers()
using LinearAlgebra
using Plots
using LaTeXStrings

#include("../../src/CTDiffFlow.jl")
#using .CTDiffFlow
using CTDiffFlow
#
# Definition of the second member
include("../../../test/ode_examples.jl")
t0 = bruss.t0
tf = bruss.tf
x0 = bruss.x0
λ = bruss.λ
rhs_bruss = bruss.rhs

tol = 1.e-4

plt1 = plot();
plt2 = plot()

Tol = 1.e-3
reltol = Tol
abstol = Tol
#algo = dopri5()
algo = Tsit5()
lw = 2

function graph_end!(plt1, plt2, δh)
    ∂x0_flow_end = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ;var_ind = :end, wrt = :x0,δh)
    ∂λ_flow_end = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ;var_ind = :end, wrt = :λ,δh)

    #println("∂x0_flow_var = ", ∂x0_flow_var(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
    #println("∂λ_flow_var = ", ∂λ_flow_var(t0, x0, tf, λ; reltol=reltol, abstol=abstol))
    Λ = range(2.88; stop=3.08, length=1001)
    n = 2;
    N = length(Λ)
    fdiff = zeros(N, n)
    ∂_λ_x0 = [0.,1]
    x0₁ = x0[1]
    for i in 1:N
        x0 = [x0₁, Λ[i]]
        fdiff[i, :] = ∂λ_flow_end(t0, x0, tf, [Λ[i]]; alg=algo, reltol=reltol, abstol=abstol) +
                      ∂x0_flow_end(t0, x0, tf, [Λ[i]]; alg=algo, reltol=reltol, abstol=abstol)*∂_λ_x0
    end
    plot!(
        plt1,
        Λ,
        fdiff[:, 1];
        xlabel="λ",
        ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)",
        lw=lw,
    )
    plot!(
        plt2,
        Λ,
        fdiff[:, 2];
        xlabel="λ",
        ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)",
        lw=lw,
    )

    #plt = plot!(plt,p1,p2,layout=(2,2),legend=false)

    #| label: fig-finite-diff-DP5
    #| fig-cap: "Derivative computing by finite differences. $t_f=20, \\lambda$ ranging from 2.88 to 3.08, $Tol=RelTol=AbsTol=10^{-4}$. Top graphs is for  $\\delta\\lambda=4Tol$ and bottom graphs for $\\delta\\lambda=\\sqrt{Tol}$. The numerical integrattion is done with DP5()."
    #algo = DP5()

end

graph_end!(plt1, plt2,4*Tol)
plt3 = plot(); plt4 = plot();
graph_end!(plt3, plt4,:default)
plt = plot(plt1, plt2, plt3, plt4)

savefig(plt, "docs/src/assets/plot_end_bruss.png")
