function test_CTFlows()
    # Tests of derivative of CTFlows
    # First example
    @testset "linear example" begin
          backend = AutoForwardDiff()
          reltol = 1.e-8
          abstol = 1.e-12
          tol_error = eps()
          atol = 4*eps()
          rtol = 4*eps()
          t0 = ode_lin1.t0
          tf = ode_lin1.tf
          x0 = ode_lin1.x0
          λ  = ode_lin1.λ
          rhs_lin1 = ode_lin1.rhs
          # end
          # wrt x0
          wrt = :x0
          #sol_∂x0_flow = ode_lin1.sol_∂x0_flow
          
          vf1 = CTFlows.VectorField(
            (x,λ) -> A(λ)*x;
            is_autonomous=true,
            is_variable=true
            )
           # IND of CTFlow
           flow_vf1 = CTFlows.Flow(vf1;reltol=reltol, abstol=abstol)
           ind_CTFlows_vf1 = jacobian(x0 -> flow_vf1(t0,x0,tf;variable=λ), backend, x0)

           # IND from CFDiffFlow
           ∂x0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
           sol_∂x0_ind = ∂x0_flow(t0, x0, tf, λ;reltol=reltol, abstol=abstol)
           @test isapprox(ind_CTFlows_vf1, sol_∂x0_ind, atol=atol, rtol=rtol)

        # Variational equation and IND from CFDiffFlow
           ∂x0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
           sol_∂x0_var = ∂x0_flow(t0, x0, tf, λ;reltol=reltol, abstol=abstol)
           @test isapprox(ind_CTFlows_vf1, sol_∂x0_var, atol=atol, rtol=rtol)
        #end;
end

end

