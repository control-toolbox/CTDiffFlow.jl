using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings
using BenchmarkTools
using ForwardDiff: ForwardDiff
using Zygote: Zygote
using SciMLSensitivity
using DifferentiationInterface

include("../../src/DiffFlow.jl")
using .DiffFlow
#
# Definition of the second member
function bruss(x, par, t)
    λ = par[1]
    x₁ = x[1]
    x₂ = x[2]
    return [1+x₁^2*x₂-(λ+1)*x₁, λ*x₁-x₁^2*x₂]
end

t0 = 0.0;
tf = 20.0
tspan = (t0, tf)
λ = [3.0]
x01 = 1.3
funx0(λ) = [x01, λ[1]]
tol = 1.e-4

plt1 = plot();
plt2 = plot()

reltol = 1.e-4;
abstol = 1.e-4

function graph_eq_var(plt1, plt2)
    ∂λ_flow_var = DiffFlow.build_∂λ_flow_var(bruss, t0, funx0, tf, λ)

    println(∂λ_flow_var(t0, funx0, tf, λ; reltol=reltol, abstol=abstol))
    # ne fonctionne pas 
    # @btime ∂λ_flow_var(t0,funx0,tf,λ;reltol=reltol,abstol=abstol)
    #println(∂x0_flow_var(t0,funx0,tf,λ;reltol=reltol, abstol=abstol))

    Λ = range(2.88; stop=3.08, length=1001)
    n = 2;
    N = length(Λ)
    fdiff = zeros(N, n)

    for i in 1:N
        fdiff[i, :] = ∂λ_flow_var(t0, funx0, tf, [Λ[i]]; reltol=reltol, abstol=abstol)
    end
    plot!(
        plt1,
        Λ,
        fdiff[:, 1];
        xlabel="λ",
        ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plot!(
        plt2,
        Λ,
        fdiff[:, 2];
        xlabel="λ",
        ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )

    #plt = plot!(plt,p1,p2,layout=(2,2),legend=false)

    #| label: fig-finite-diff-DP5
    #| fig-cap: "Derivative computing by finite differences. $t_f=20, \\lambda$ ranging from 2.88 to 3.08, $Tol=RelTol=AbsTol=10^{-4}$. Top graphs is for  $\\delta\\lambda=4Tol$ and bottom graphs for $\\delta\\lambda=\\sqrt{Tol}$. The numerical integrattion is done with DP5()."
    #algo = DP5()

end

graph_eq_var(plt1, plt2)
plt = plot(plt1, plt2)
savefig(plt, "article/figures/plot_var1.png")

function graph_diff_auto_flow(plt1, plt2)
    ∂λ_flow = DiffFlow.build_∂λ_flow(bruss, t0, funx0, tf, λ)
    δλ = [1.0]
    println(∂λ_flow(t0, funx0, tf, λ; reltol=reltol, abstol=abstol))

    Λ = range(2.88; stop=3.08, length=1001)
    n = 2;
    N = length(Λ)
    fdiff = zeros(N, n)

    for i in 1:N
        fdiff[i, :] = ∂λ_flow(t0, funx0, tf, [Λ[i]]; reltol=reltol, abstol=abstol)
    end
    plot!(
        plt1,
        Λ,
        fdiff[:, 1];
        color=:magenta,
        xlabel="λ",
        ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plot!(
        plt2,
        Λ,
        fdiff[:, 2];
        color=:magenta,
        xlabel="λ",
        ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )

    #plt_diff_flow1 = plot(p1,p2,layout=(2,2),color=:magenta,legend=false)
    #plt = plot!(plt,p1,p2,layout=(2,2),color=:magenta,legend=false)

    #| label: fig-finite-diff-DP5
    #| fig-cap: "Derivative computing by finite differences. $t_f=20, \\lambda$ ranging from 2.88 to 3.08, $Tol=RelTol=AbsTol=10^{-4}$. Top graphs is for  $\\delta\\lambda=4Tol$ and bottom graphs for $\\delta\\lambda=\\sqrt{Tol}$. The numerical integrattion is done with DP5()."
    algo = DP5()

    #savefig(plt_diff_flow1, "plot_diff_flow1.png")
end

graph_diff_auto_flow(plt1, plt2)

plt = plot(plt1, plt2)
savefig(plt, "article/figures/plot_diff_flow1.png")
