using OrdinaryDiffEq,  Test, DiffEqDevTools,
      ODEInterface, ODEInterfaceDiffEq

import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

prob = prob_ode_linear
#=
prob.tspan
prob.f
prob.u0
prob.p
prob.kwargs
prob.problem_type
=#

sol = solve(prob, DP5())

sol2 = solve(prob, dopri5())
println("sol.t[2] - sol2.t[2] = ", sol.t[2] - sol2.t[2])
@test sol.t[2] ≈ sol2.t[2]

prob = prob_ode_2Dlinear
sol = solve(prob, DP5(), internalnorm = (u, t) -> sqrt(sum(abs2, u)))

# Change the norm due to error in dopri5.f
sol2 = solve(prob, dopri5())
println("sol.t[2] - sol2.t[2] = ", sol.t[2] - sol2.t[2])
@test sol.t[2] ≈ sol2.t[2]

prob = deepcopy(prob_ode_linear)
prob
#prob2 = ODEProblem(prob.f, prob.u0, (1.0, 0.0), 1.01)
prob2 = ODEProblem(prob.f, prob.u0, (0.0, 1.0), 1.01)
sol = solve(prob2, DP5())

sol2 = solve(prob2, dopri5())
println("sol.t[2] - sol2.t[2] = ", sol.t[2] - sol2.t[2])
@test sol.t[2] ≈ sol2.t[2]

prob = deepcopy(prob_ode_2Dlinear)
#prob2 = ODEProblem(prob.f, prob.u0, (1.0, 0.0), 1.01)
prob2 = ODEProblem(prob.f, prob.u0, (0.0, 1.0), 1.01)
sol = solve(prob2, DP5(), internalnorm = (u, t) -> norm(u)/sqrt(length(u)))#sqrt(sum(abs2, u)))

# Change the norm due to error in dopri5.f
sol2 = solve(prob2, dopri5())
println("sol.t[2] - sol2.t[2] = ", sol.t[2] - sol2.t[2])
@test sol.t[2] ≈ sol2.t[2]