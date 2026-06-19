function test_end_ind_var()
# unitary test for the derivative of the flow
  tol_error = eps()
  @testset verbose = true "end_ind_var" begin
        backend = AutoForwardDiff()
        # linear ode example
        @testset "linear example" begin
          reltol = 1.e-8
          abstol = 1.e-12
          t0 = ode_lin1.t0
          tf = ode_lin1.tf
          x0 = ode_lin1.x0
          λ  = ode_lin1.λ
          rhs_lin1 = ode_lin1.rhs
          # end
          # wrt x0
          wrt = :x0
          sol_∂x0_flow = ode_lin1.sol_∂x0_flow
          atol = sqrt(abstol)
          rtol = sqrt(reltol)
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂x0_flow, sol_end, atol=atol, rtol=rtol)
          # ind
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂x0_flow, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂x0_flow, sol_var, atol=atol, rtol=rtol)   
          atol = eps()
          rtol = eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)   
          # wrt λ
          wrt = :λ
          sol_∂λ_flow = ode_lin1.sol_∂λ_flow
          atol = sqrt(abstol)
          rtol = sqrt(reltol)
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂λ_flow, sol_end, atol=atol, rtol=rtol)
          # ind
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂λ_flow, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂λ_flow, sol_var, atol=atol, rtol=rtol)   
          atol = eps()
          atol = eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol) 
          # wrt t0
          wrt = :t0
          sol_∂λ_flow = ode_lin1.sol_∂t0_flow
          atol = sqrt(abstol)
          rtol = sqrt(reltol)
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂t0_flow, sol_end, atol=atol, rtol=rtol)
          # ind
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂t0_flow, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂t0_flow, sol_var, atol=atol, rtol=rtol)   
          atol = eps()
          atol = eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)   
          # wrt tf
          wrt = :tf
          sol_∂λ_flow = ode_lin1.sol_∂tf_flow
          atol = sqrt(abstol)
          rtol = sqrt(reltol)
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂tf_flow, sol_end, atol=atol, rtol=rtol)
          # ind
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂tf_flow, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_lin1, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_∂tf_flow, sol_var, atol=atol, rtol=rtol)   
          atol = eps()
          atol = eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)  
        end
        
        # brusselator example
        @testset "brusselator example" begin
          t0 = bruss.t0
          tf = bruss.tf
          x0 = bruss.x0
          λ = bruss.λ
          rhs_bruss = bruss.rhs
          reltol = 1.e-8
          abstol = 1.e-8
          # wrt x0
          wrt = :x0
          rtol = sqrt(reltol)
          atol = sqrt(abstol)
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          # ind
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_end, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂x0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂x0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          atol = 10*eps()
          rtol = 10*eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)      
          # wrt λ
          wrt = :λ
          rtol = 4*sqrt(reltol)
          atol = 4*sqrt(abstol)
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          # ind
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_end, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂λ_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂λ_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          atol = 10*eps()
          rtol = 10*eps()
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)    
          # wrt t0
          wrt = :t0
          rtol = 4*sqrt(reltol)
          atol = 4*sqrt(abstol)
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          # ind
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_end, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂t0_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂t0_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          #atol = 10*eps()
          #rtol = 10*eps()
          @test isapprox(sol_end, sol_var, atol=atol, rtol=rtol)    
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)  
          # wrt tf
          wrt = :tf
          rtol = 4*sqrt(reltol)
          atol = 4*sqrt(abstol)
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :end, wrt = wrt)#, δh=1.e-10)
          sol_end = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          # ind
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :ind, wrt = wrt)
          sol_ind = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          @test isapprox(sol_end, sol_ind, atol=atol, rtol=rtol)      
          # var
          ∂tf_flow = CTDiffFlow.build_∂flow(rhs_bruss, t0, x0, tf, λ; var_ind= :var, wrt = wrt)
          sol_var = ∂tf_flow(t0, x0, tf, λ; adaptive=true, reltol=reltol, abstol=abstol)
          #atol = 10*eps()
          #rtol = 10*eps()
          @test isapprox(sol_end, sol_var, atol=atol, rtol=rtol)    
          @test isapprox(sol_ind, sol_var, atol=atol, rtol=rtol)  
        end
        
      end;
    end