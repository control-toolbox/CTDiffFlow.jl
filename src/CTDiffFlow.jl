module CTDiffFlow

using DifferentiationInterface
using OrdinaryDiffEq
using LinearAlgebra

include("./myode43/myode43.jl")

struct Sol
    t
    u
end
# Derivative with respect to x0
function build_∂x0_flow(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; var_ind = :var, backend = AutoForwardDiff()) #,print_step=false)
    """
        Return the function which compute the derivative of the flow
        
        input
        -----
        rhs : right hand side which defines the flow
              xpoint = rhs(x,λ,t)
        t0 : the initial time
             Real
        x0 : initial state
             Vector(Real)
        tf : final time
             Real
        λ : parameter
            Vector of Real
        var_ind : Integrate the variational equations or Internal Numerical Equation (Automatic differentiation on the flow)
              var_ind = :var -> Variational equations
              vat_ind = :ind -> IND
        backend : backend of the DifferentiationInterface.jl

        output
        ------
        the function ∂x0_flow the derivative with respect to ine initial condition of the flow 
        at the final time tf or on the interval (t0,tf)
        ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        or
        ∂x0_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
    """
    @assert (var_ind == :var || var_ind == :ind) "Error the var_ind optional argument of the build_∂x0_flow function is not equal to :var or :ind"
    if var_ind == :var
        build_∂x0_flow_var(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = backend)
    elseif var_ind == :ind
        build_∂x0_flow_ind(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = backend)
    end
end

function built_rhs_var(rhs::Function; wrt = :x0, backend = AutoForwardDiff() )
    """
        Built the second member of the variational equations
        input
        -----
        rhs : right hand side
              xpoint = rhs(x,λ,t)
        wrt : with respect to
              wrt = :x0 -> with respect to the intial condition
              wrt = :λ -> with respect to the parameter λ
              wrt = :t0 -> with respect to the initial time
        backend : backend of the DifferentiationInterface.jl

        output
        ------
        the fucntion rhs_var(xX,λ,t)
            xX is a matrix (n,n+1) for wrt=:x0, (n,p) for wrt=:λ, (n,1) for wrt=:t0
            The fisrt column is the initial state variable x
            The 2:end column is the state X of the variational equations
    """

    @assert (wrt == :x0 || wrt == :λ || wrt ==:t0) "Error the wrt optional argument of the built_rhs_var function is not equal to :x0 or :λ ot :t0"
    fun_x(x,λ,t) = jacobian(x -> rhs(x,λ,t), backend, x)
    if wrt == :λ
        fun_λ(x,λ,t) = jacobian(λ -> rhs(x,λ,t), backend, λ)
        function rhs_var_λ(xX,λ,t)
            x = xX[:,1]
            X = xX[:,2:end]
            xpoint = rhs(x,λ,t)
            Xpoint = fun_x(x,λ,t)*X + fun_λ(x,λ,t)
            return [xpoint Xpoint]
        end
        return rhs_var_λ
    else
        function rhs_var_x0(xX,λ,t)
            x = xX[:,1]
            X = xX[:,2:end]
            xpoint = rhs(x,λ,t)
            Xpoint = fun_x(x,λ,t)*X
            return [xpoint Xpoint]
        end
           return rhs_var_x0
    end
end

# Integration of the variational equations
# derivatives with respect to x0
function build_∂x0_flow_var(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff()) #,print_step=false)

    rhs_var = built_rhs_var(rhs , wrt = :x0, backend = backend)


    function ∂x0_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        n = length(x0)
        x0δx0 = [x0 Matrix(I(n))]
        algo = get(ode_kwargs, :alg, Tsit5())
        reltol = get(ode_kwargs, :reltol, 1.e-3)
        abstol = get(ode_kwargs, :abstol, 1.e-6)
        @assert (typeof(reltol)<:Real || length(reltol)==n || size(reltol)==(n,n+1)) "Error in the dimension of reltol" 
        adaptive = get(ode_kwargs, :adaptive, true)
        if adaptive
            my_Inf = prevfloat(typemax(Float64))
            n = length(x0)
            p = n
            if typeof(reltol) <: Real
              RelTol = [reltol*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
            elseif length(reltol)==n
              RelTol = [reltol.*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
            else  # reltol is a Matix (n,n+1)
              RelTol = reltol
            end
            if typeof(abstol) <: Real
              AbsTol = [abstol*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
            elseif length(reltol)==n
              AbsTol = [abstol.*ones(n,1) my_Inf*ones(n,p)]/sqrt(p+1)
            else  # reltol is a Matix (n,n+1)
              AbsTol = abstol
            end
            algo = get(ode_kwargs, :alg, Tsit5())
            if algo=="myode43"
              internalnorm = get(ode_kwargs, :internalnorm, :default)
              T,X = myode43(rhs_var,x0δx0,λ,(t0,tf),RelTol,AbsTol; internalnorm=internalnorm)
              sol = Sol(T,X) # for having the same structure as standard numerical integration
             return sol
            else
              ivp = ODEProblem(rhs_var, x0δx0, tspan, λ)
              sol = solve(ivp, alg=algo; ode_kwargs..., reltol=RelTol, abstol=AbsTol)
              return sol
            end
        else
            t0 = tspan[1]; tf = tspan[2];
            dt = get(ode_kwargs, :dt, (tf-t0)/100)
            if algo=="myode43"
              T,X = myode43(rhs_var,x0δx0,λ,t0:dt:tf)
              sol = Sol(T,X) # for having the same structure as standard numerical integration
             return sol
            else
              ivp = ODEProblem(rhs_var, x0δx0, tspan, λ)
              sol = solve(ivp, alg=algo; ode_kwargs...)
              return sol
            end

        end
    end

    function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        sol = ∂x0_flow((t0,tf),x0,λ; print_times, ode_kwargs...)
        if print_times
            return sol.u[end][:,2:end], sol.t
        else
            return sol.u[end][:,2:end]
        end
     end
    return ∂x0_flow
end

# derivatives with respect to λ
function build_∂λ_flow_var(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff()) #,print_step=false)
    rhs_var = built_rhs_var(rhs , wrt = :λ, backend = backend)

    function ∂λ_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, λ::Vector{<:Real}; ode_kwargs...)
        n = length(x0); p = length(λ)
        x0δλ0 = [x0 zeros(n,p)]
        ivp = ODEProblem(rhs_var, x0δλ0, tspan, λ)
        algo = get(ode_kwargs, :alg, Tsit5())
        sol = solve(ivp, alg=algo; ode_kwargs...)
        return sol
    end

    function ∂λ_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; ode_kwargs...)
        sol = ∂λ_flow((t0,tf),x0,λ; ode_kwargs...)
        return sol.u[end][:,2:end]
     end
    return ∂λ_flow
end

function build_∂x0_flow_ind(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff())

    # Jacobian matrix
    function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        
# Automatic differentiation on the flow
# derivatives with respect to x0
        algo = get(ode_kwargs, :alg, Tsit5())
        adaptive = get(ode_kwargs, :adaptive, true)
        if algo=="myode43"
            Rtol = get(ode_kwargs, :reltol, 1.e-3)
            Atol = get(ode_kwargs, :abstol, 1.e-6)
            internalnorm = get(ode_kwargs, :internalnorm, :default)
            T,X = myode43(rhs,x0,λ,(t0,tf),Rtol,Atol; internalnorm=internalnorm)
        else
            ivp = ODEProblem(rhs, x0, (t0,tf), λ)
            sol = solve(ivp, alg=algo; ode_kwargs...)
            T = sol.t
        end


        
        function _flow(x0)
            if algo=="myode43"
                if adaptive
                  Rtol = get(ode_kwargs, :reltol, 1.e-3)
                  Atol = get(ode_kwargs, :abstol, 1.e-6)
                  T,X = myode43(rhs,x0,λ,(t0,tf),Rtol,Atol)
                else
                  dt = get(ode_kwargs, :dt, (tf-t0)/100)
                  T,X = myode43(rhs,x0,λ,t0:dt:tf)
                end    
                return X[end]
            else

                ivp = ODEProblem(rhs, x0, (t0,tf), λ)
                sol = solve(ivp, alg=algo; ode_kwargs...)
                T = sol.t
                return sol.u[end]
            end
        end
        if print_times
            return jacobian(_flow,backend,x0), T
        else
            return jacobian(_flow,backend,x0)
        end
     end
    
    return ∂x0_flow
end

# derivatives with respect to λ
function build_∂λ_flow_ind(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff())
    
    function ∂λ_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        T = []
        function _flow(λ)
            ivp = ODEProblem(rhs, x0, tspan, λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            sol = solve(ivp, alg=algo; ode_kwargs...)
            T = sol.t
            
            return reduce(hcat,sol.u)
        end
        if print_times
            return jacobian(_flow,backend,λ), T
        else
            return jacobian(_flow,backend,λ)
        end
     end
     
    function ∂λ_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ0::Vector{<:Real}; print_times=false, ode_kwargs...)
        T = []
        function _flow(λ)
            ivp = ODEProblem(rhs, x0, (t0,tf), λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            sol = solve(ivp, alg=algo; ode_kwargs...)
            T = sol.t
            return sol.u[end]
        end
        if print_times
            return jacobian(_flow,backend,λ0), T
        else
            return jacobian(_flow,backend,λ0)
        end
     end
  return ∂λ_flow
end

end # CTDiffFlow
