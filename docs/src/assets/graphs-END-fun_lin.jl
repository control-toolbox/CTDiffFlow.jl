using Pkg
Pkg.activate(".")
using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings
using Test
# Definition of the second member

include("../../test/fun_examples.jl")


function graph_END(ind)
    # ind = indice of the derivative in λ
    function END(fun, fun_x0, tf, tol, λ, δλ, ind; algo=Tsit5())
        tspan = (0.0, tf)
        n = length(fun_x0(λ))
        p = length(λ)
        ∂λ_xtf = zeros(n,p)
        ivp = ODEProblem(fun, fun_x0(λ), tspan, (λ))
        sol1 = solve(ivp; alg=algo, reltol=tol, abstol=tol)
            ej = zeros(p); ej[ind] = 1;
            x0 = funx0(λ+δλ*ej)
            ivp = ODEProblem(fun, x0, tspan, λ+δλ*ej)
            sol2 = solve(ivp; alg=algo, reltol=tol, abstol=tol)
            ∂λ_xtf = (sol2.u[end]-sol1.u[end])/δλ
        return ∂λ_xtf
    end


    λ = [1.,2]; p = length(λ);
    t0 = 0. ; tf = 1.; tspan = (t0,tf);
    funx0(λ) = [λ[2],1.,1]; 
    n = length(funx0(λ))
    tol = 1.e-10
    δλ = sqrt(tol)
    # test of ∂λ_xtf
    tf_test = 1.
    @test isapprox(sol_∂λ_flow(tf_test,λ)[:,1], END(fun_lin, fun_lin_x0, tf_test, tol, λ, δλ,1; algo=Tsit5()), atol = 10*δλ)


    Λ = range(0.5*λ[ind]; stop=1.5*λ[ind], length=11)

    N = length(Λ)
    fdiff = zeros(N, n)
    ∂λ_flow = zeros(N, n)
    tol = 1.e-4
    δλ = tol
    for i in 1:N
        λi = λ
        λi[ind] = Λ[i]
        fdiff[i, :] = END(fun_lin, fun_lin_x0, tf, tol, λi, δλ,ind)
        ∂λ_flow[i,:] = sol_∂λ_flow(tf,λi)[:,ind]
    end
    p1 = plot(
        Λ,
        fdiff[:, 1];
        xlabel="λ",
        ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plot!(p1,Λ,∂λ_flow[:,1])
    p2 = plot(
        Λ,
        fdiff[:, 2];
        xlabel="λ",
        ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plot!(p2,Λ,∂λ_flow[:,2])

      p3 = plot(
        Λ,
        fdiff[:, 3];
        xlabel="λ",
        ylabel=L"\frac{\partial x_3}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plot!(p3,Λ,∂λ_flow[:,3])

    δλ = sqrt(tol)
    for i in 1:N
        λi = λ
        λi[ind] = Λ[i]
        fdiff[i, :] = END(fun_lin, fun_lin_x0, tf, tol, λi, δλ, ind)
    end
    p4 = plot(
        Λ,
        fdiff[:, 1];
        xlabel="λ",
        ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    p5 = plot(
        Λ,
        fdiff[:, 2];
        xlabel="λ",
        ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
         p6 = plot(
        Λ,
        fdiff[:, 3];
        xlabel="λ",
        ylabel=L"\frac{\partial x_3}{\partial \lambda}(t_f,\lambda)",
        lw=3,
    )
    plt_END = plot(p1, p2, p3, p4, p5, p6; layout=(3, 2), legend=false)


    #plt_END = plot(p5, p6, p7, p8; layout=(2, 2), legend=false)


    return plt_END
end

plt_END = graph_END(1)
savefig(plt_END, "docs/assets/plot_END_1.png")

plt_END = graph_END(2)
savefig(plt_END, "docs/assets/plot_END_2.png")
