###############################################################################
# From https://github.com/SciML/ODE.jl/blob/master/src/ODE.jl
# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
#
#function hinit(F, x0, t0::T, tend, p, reltol, abstol) where T
function hinit(F, x0, t0, tend, p, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend-t0)
    tdir==0 && error("Zero time span")
    tau = max(reltol*norm(x0, Inf), abstol)
    d0 = norm(x0, Inf)/tau
    f0 = F(t0, x0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = x0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max((10)^(-6), (10)^(-3)*h0)
    else
        pow = -(2 + log10(max(d1, d2)))/(p + 1)
        h1 = 10^pow
    end

    return tdir*min(100*h0, h1, tdir*(tend-t0)), tdir, f0
end

A(λ) = [λ[1] 0 0; 0 λ[2] 0; 0 0 λ[1]-λ[2]]
fun1(t, x) = A(λ)*x

function fun2(t, xX)
    [fun1(xX[:, 1], λ, t) fun1(
        xX[:, 2:end], λ, t
    )+[xX[1, 1] 0; 0 xX[2, 1]; xX[3, 1] xX[3, 1]]]
end

t0 = 0.0;
tend = 2.0;
p = 2
x0 = [0.0, 1.0, 1]
reltol = 1.e-4;
abstol = 1.e-4

hinit(fun1, x0, t0, tend, p, reltol, abstol)

xX0 = [funx0(λ) [0 1; 0 0; 0 0]]
hinit(fun2, xX0, t0, tend, p, reltol, abstol)
RelTol = reltol*ones(n, p+1)
AbsTol = RelTol

hinit(fun2, xX0, t0, tend, p, RelTol, AbsTol)

RelTol = [reltol*ones(n, 1) Inf*ones(n, p)]/sqrt(p+1)
AbsTol = [abstol*ones(n, 1) Inf*ones(n, p)]/sqrt(p+1)
