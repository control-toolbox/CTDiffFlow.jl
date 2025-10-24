using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings
# Definition of the second member
function bruss(x,par,t)
  λ = par[1]
  x₁ = x[1]
  x₂ = x[2]
  return [1+x₁^2*x₂-(λ+1)*x₁ , λ*x₁-x₁^2*x₂]
end

tf = 20
tspan = (0.,tf)
x₁0 = 1.3
par = [3.]
tol = 1.e-4

function graph_END()
function END(fun,x₁0,tf,tol,λ,δλ;algo=Tsit5())
reltol = tol; abstol = tol;
  tspan = (0.,tf)
  x0 = [x₁0,λ+δλ]
  ivp = ODEProblem(fun, x0, tspan, (λ+δλ))
  sol2 = solve(ivp, alg=algo, reltol = tol, abstol = tol)
  x0 = [x₁0,λ]
  ivp = ODEProblem(fun, x0, tspan, (λ))
  sol1 = solve(ivp, alg=algo, reltol = tol, abstol = tol)
  return (sol2.u[end]-sol1.u[end])/δλ
end

Λ = range(2.88,stop=3.08,length=1001)
n = 2; N = length(Λ)
fdiff = zeros(N,n)
δλ = 4*tol
for i in 1:N
  fdiff[i,:] = END(bruss,x₁0,tf,tol,Λ[i],δλ)
end
p1 = plot(Λ,fdiff[:,1],xlabel="λ", ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)", lw = 3)
p2 = plot(Λ,fdiff[:,2],xlabel="λ", ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)", lw = 3)

δλ = sqrt(tol)
for i in 1:N
  fdiff[i,:] =END(bruss,x₁0,tf,tol,Λ[i],δλ)
end
p3 = plot(Λ,fdiff[:,1],xlabel="λ", ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)", lw = 3)
p4 = plot(Λ,fdiff[:,2],xlabel="λ", ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)", lw = 3)
plot(p1,p2,p3,p4,layout=(2,2),legend=false)

#| label: fig-finite-diff-DP5
#| fig-cap: "Derivative computing by finite differences. $t_f=20, \\lambda$ ranging from 2.88 to 3.08, $Tol=RelTol=AbsTol=10^{-4}$. Top graphs is for  $\\delta\\lambda=4Tol$ and bottom graphs for $\\delta\\lambda=\\sqrt{Tol}$. The numerical integrattion is done with DP5()."
#algo = DP5()
algo = Tsit5()
δλ = 4*tol
for i in 1:N
  fdiff[i,:] = END(bruss,x₁0,tf,tol,Λ[i],δλ;algo=algo)
end
p5 = plot(Λ,fdiff[:,1],xlabel="λ", ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)", lw = 3)
p6 = plot(Λ,fdiff[:,2],xlabel="λ", ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)", lw = 3)

δλ = sqrt(tol)
algo = DP5()
for i in 1:N
  fdiff[i,:] = END(bruss,x₁0,tf,tol,Λ[i],δλ;algo=algo)
end
p7 = plot(Λ,fdiff[:,1],xlabel="λ", ylabel=L"\frac{\partial x_1}{\partial \lambda}(t_f,\lambda)", lw = 3)
p8 = plot(Λ,fdiff[:,2],xlabel="λ", ylabel=L"\frac{\partial x_2}{\partial \lambda}(t_f,\lambda)", lw = 3)
plt_END = plot(p5,p6,p7,p8,layout=(2,2),legend=false)

savefig(plt_END, "article/figures/plot_END.png")
end

graph_END()

