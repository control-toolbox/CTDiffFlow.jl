function test_eq_var_rhs()
# unitary test for the derivative of the right hand side of the ivp 
# wrt the initial condition x0 and the parameter λ
  tol_error = eps()
  @testset verbose = true "eq_var_rhs" begin
        backend = AutoForwardDiff()
        # linear ode example
        @testset "linear example" begin
          t0 = ode_lin1.t0
          x0 = ode_lin1.x0
          λ = ode_lin1.λ
          xδx = [x0 I]
          rhs_var_x0 =  CTDiffFlow.built_rhs_var(ode_lin1.rhs; wrt = :x0, backend = backend )
          @test isapprox(rhs_var_x0(xδx,λ,t0), ode_lin1.eq_var_x0(xδx,λ,t0), atol=tol_error)
          rhs_var_λ =  CTDiffFlow.built_rhs_var(ode_lin1.rhs; wrt = :λ, backend = backend )
          xδλ = xδx[:,1:3]
          @test isapprox(rhs_var_λ(xδλ,λ,t0), ode_lin1.eq_var_λ(xδλ,λ,t0), atol=tol_error)
        end
        #
        # brusselator example
        @testset "brusselator example" begin
          t0 = bruss.t0
          x0 = bruss.x0
          λ = bruss.λ
          xδx = [x0 I]
          rhs_var_x0 =  CTDiffFlow.built_rhs_var(bruss.rhs; wrt = :x0, backend = backend )
          @test isapprox(rhs_var_x0(xδx,λ,t0), bruss.eq_var_x0(xδx,λ,t0), atol=tol_error)
          rhs_var_λ =  CTDiffFlow.built_rhs_var(bruss.rhs; wrt = :λ, backend = backend )
          p = length(λ)
          xδλ = xδx[:,1:p+1]
          @test isapprox(rhs_var_λ(xδλ,λ,t0), bruss.eq_var_λ(xδλ,λ,t0), atol=tol_error)
        end
      end;
    end