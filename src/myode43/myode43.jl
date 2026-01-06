#
# Numerical integration with control step size
# algorithm : see Hairer Tome 1 page
#
#
using LinearAlgebra
function inith(rhs, x0, par, t0, Atol, Rtol)
    #
    #
    n = length(x0);
    sc = Atol .+ abs.(x0) .* Rtol;
    println("sc = ", sc)
    println("x0 = ", x0)
    k0 = rhs(x0, par, t0);
    println(x0 ./ sc)
    println(norm(x0 ./ sc))
    d0 = norm(x0 ./ sc)/sqrt(n);              # normalement c'est la bonne valeur cf page 168
    d1 = norm(k0 ./ sc)/sqrt(n);              # normalement c'est la bonne valeur cf page 168
    #d0 =  norm(x0./sc); 
    #d1 =  norm(k0./sc);
    h0 = 0.01*(d0/d1);
    if (d0 < 1.e-5) || (d1 < 1.e-5)
        h0 = 1.e-6;
    end;
    x1 = x0 + h0*k0;
    k1 = rhs(x1, par, t0+h0);
    d2 = norm((k1-k0) ./ sc)/sqrt(n)/h0;      # normalement c'est la bonne valeur cf page 168
    #d2 =  norm((k1-k0)./sc)/h0;
    if max(d1, d2) < 1.e-15
        h1 = max(1.e-6, h0*1.e-3);
    else
        h1 = (0.01/max(d1, d2))^(1/4);
    end
    h = min(100*h0, h1);
    return h
end

function myode43(rhs, x0, par, t0tf, Rtol, Atol)
    #
    # Initialisation
    Scal_Type=eltype(par)
    p=4; # ordre
    t0=t0tf[1];
    tf=t0tf[2];
    hmax=tf-t0;
    Npasmax=1000;
    n = length(x0);
    T = [Scal_Type(t0)];
    #T = [t0];
    X = [Scal_Type.(x0)]
    #
    #
    # step initialisation
    #h=0.03;
    h=inith(rhs, x0, par, t0, Atol, Rtol)
    #println("h= ",h)
    #
    nstep=0;
    t=t0
    x=x0;
    fin=0;
    while (nstep < Npasmax) && (fin == 0)
        nstep = nstep+1
        # Runge-Kutta method
        k1 = rhs(x, par, t);
        k2 = rhs(x+(h/3)*k1, par, t+h/3);
        k3 = rhs(x+h*(-k1/3+k2), par, t+2*h/3);
        k4 = rhs(x+h*(k1-k2+k3), par, t+h);
        x1 = x+(h/8)*(k1+3*k2+3*k3+k4);
        xhat1 = x+(h/12)*(k1+6*k2+3*k3+2*rhs(x1, par, t));
        # err
        sc = Atol .+ max.(abs.(x), abs.(x1)) .* Rtol;
        #println("sc = ", sc)
        err = norm((x1-xhat1) ./ sc)/sqrt(n);
        # println("err = ", err)
        #
        # calcul du pas
        if (err < 1)
            t = t+h;
            x = x1;
            #println("t= ",t)
            #println("T= ",T)
            #println("X= ",X)
            push!(T, t)
            #println("T= ",T)
            #println("x1= ",x1)
            push!(X, x1)
            #println("X= ",X)
            if (t > tf-h/2)
                fin = 1;
            end;
        end;
        h = h*min(5, max(0.2, 0.9*(1/err)^(1/p)));
        if (t+h > tf)
            h = tf-t;
        end;
    end;

    return T, X
end

function get_Xij(X, i, j)
    """ 
        Get from the vector of vector X the i,j composante for all times
    """
    N = length(X)
    Xij = zeros(N)
    for l in 1:N
        Xij[l] = X[l][i, j]
    end
    return Xij
end
