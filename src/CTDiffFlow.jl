module CTDiffFlow

using DifferentiationInterface
using OrdinaryDiffEq
using LinearAlgebra

function built_rhs_var(rhs::Function; wrt = :x0, backend = AutoForwardDiff() )
    """
        Built the second member of the variational equations
        input
        -----
        rhs : right hand side
            xpoint = rhs(x,λ,t)
        wrt : with respect to
              wrt = :x0 -> with respevt to the intial condition
              wrt = :λ -> with respect to the parameter λ
              wrt = :t0 -> with respect to the initial time
        backend : backend of the DifferentiationInterface.jl

        output
        ------
        the fucntion rhs_var(xX,λ,t)
            xX is a matrix (n,n+1) 
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
        ivp = ODEProblem(rhs_var, x0δx0, tspan, λ)
        algo = get(ode_kwargs, :alg, Tsit5())
        sol = solve(ivp, alg=algo; ode_kwargs...)
        return sol
    end

    function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; ode_kwargs...)
        #println("reltol = ", get(ode_kwargs, :reltol, 1.e-3))
        #println("abstol = ", get(ode_kwargs, :abstol, 1.e-6))
        sol = ∂x0_flow((t0,tf),x0,λ; ode_kwargs...)
        return sol.u[end][:,2:end]
     end
    return ∂x0_flow
end

# derivatives with respect to λ
function build_∂λ_flow_var(rhs::Function,t0::Real,x0::Function,tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff()) #,print_step=false)
    Jacλx0 = jacobian(x0,backend,λ)
    rhs_var = built_rhs_var(rhs , wrt = :λ, backend = backend)

    function ∂λ_flow(tspan::Tuple{<:Real,<:Real},x0::Function, λ::Vector{<:Real}; ode_kwargs...)
        x0λ = x0(λ)
        #=
        n = length(x0λ); p = length(λ);
        =#
        x0δλ0 = [x0λ Jacλx0]
        ivp = ODEProblem(rhs_var, x0δλ0, tspan, λ)
        algo = get(ode_kwargs, :alg, Tsit5())
        sol = solve(ivp, alg=algo; ode_kwargs...)
        return sol
    end

    function ∂λ_flow(t0::Real,x0::Function, tf::Real, λ::Vector{<:Real}; ode_kwargs...)
        #println("reltol = ", get(ode_kwargs, :reltol, 1.e-3))
        #println("abstol = ", get(ode_kwargs, :abstol, 1.e-6))
        sol = ∂λ_flow((t0,tf),x0,λ; ode_kwargs...)
        return sol.u[end][:,2:end]
     end
    return ∂λ_flow
end

# Automatic differentiation on the flow
# derivatives with respect to x0
function build_∂x0_flow(rhs::Function,t0::Real,x0::Vector{<:Real},tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff())

    #=
    function ∂x0_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, δx0::Vector{<:Real}, λ::Vector{<:Real}; ode_kwargs...)
        function flow(x0)
            ivp = ODEProblem(rhs, x0, tspan, λ)
            sol = solve(ivp, alg=algo; ode_kwargs...)
            return sol.u
        end
        println("coucou")
        return jacobian(flow,backend,x0)
    end
    =#
#=
       function ∂x0_flow(tspan::Tuple{<:Real,<:Real},x0::Vector{<:Real}, δx0::Matrix{<:Real}, λ::Vector{<:Real}; ode_kwargs...)
        function flow(x0)
            ivp = ODEProblem(rhs, x0, tspan, λ)
            sol = solve(ivp, alg=algo; ode_kwargs...)
            return sol.u[end]
        end
        return jacobian(flow,backend,x0)*δx0
    end
    =#

    
    function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        function _flow(x0)
            ivp = ODEProblem(rhs, x0, (t0,tf), λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            sol = solve(ivp, alg=algo; ode_kwargs...)
            return sol.u[end]
        end
        if print_times
            return jacobian(_flow,backend,x0), sol.t
        else
            return jacobian(_flow,backend,x0)
        end
     end
    
     function ∂x0_flow(t0::Real,x0::Vector{<:Real}, tf::Real, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        function _flow(x0)
            ivp = ODEProblem(rhs, x0, (t0,tf), λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            sol = solve(ivp, alg=algo; ode_kwargs...)
            return sol.u[end]
        end
        if print_times
            return jacobian(_flow,backend,x0), sol.t
        else
            return jacobian(_flow,backend,x0)
        end
     end
    ∂x0_flow
end

# derivatives with respect to λ
function build_∂λ_flow(rhs::Function,t0::Real,x0::Function,tf::Real, λ::Vector{<:Real}; backend = AutoForwardDiff())
    
    function ∂λ_flow(tspan::Tuple{<:Real,<:Real},x0::Function, λ::Vector{<:Real}; print_times=false, ode_kwargs...)
        T = []
        function _flow(λ)
            ivp = ODEProblem(rhs, x0(λ), tspan, λ)
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
    
    function ∂λ_flow(t0::Real,x0::Function, tf::Real, λ0::Vector{<:Real}; print_times=false, ode_kwargs...)
        T = []
        function _flow(λ)
            ivp = ODEProblem(rhs, x0(λ), (t0,tf), λ)
            algo = get(ode_kwargs, :alg, Tsit5())
            println("algo = ", algo)
            sol = solve(ivp, alg=algo; ode_kwargs...)
            T = sol.t
            println(sol.u[end])
            return sol.u[end]
        end
        if print_times
            println("backend = ", backend)
            return jacobian(_flow,backend,λ0), T
        else
            return jacobian(_flow,backend,λ0)
        end
     end
  return ∂λ_flow
end

end # CTDiffFlow
