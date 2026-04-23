#
# Numerical integration with control step size
# algorithm : see Hairer Tome 1 page
#
#
using LinearAlgebra
using ForwardDiff: ForwardDiff 

# variable steps
function myode43(rhs,x0,par,t0tf,Rtol,Atol;internalnorm=:default)

# Definition of the norm for having the good steps
# when datas are ForwardDiff.Dual numbers
if internalnorm == :default
    sse(x::Number) = x^2
    sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) #+ sum(sse, ForwardDiff.partials(x))
    totallength(x::Number) = 1
    function totallength(x::ForwardDiff.Dual)
      totallength(ForwardDiff.value(x)) #+ sum(totallength, ForwardDiff.partials(x)).  
    end
    totallength(x::AbstractArray) = sum(totallength, x)
    function mynorm(u)
      return sqrt(sum(x -> sse(x), u) / totallength(u))
    end
else
  # internal norm as in DifferentialEquation.jl
  t = t0tf[1]
  mynorm = u -> internalnorm(u,t)
end

#=
if eltype(x0) <: ForwardDiff.Dual
  Atol = Atol*
end
=#

function inith(rhs,x0,par,t0,Atol,Rtol)
  #
  n = length(x0);
  #println("Atol = ", Atol)
  #println("Rtol = ", Rtol)
  #println("abs.(x0) .* Rtol =" , abs.(x0) .* Rtol)
  sc = Atol .+ abs.(x0) .* Rtol;
  #println("sc = ", sc)
  k0 = rhs(x0, par, t0);
  #println("k0 = ",k0)
  #println("x0./sc = ", x0./sc)
  d0 =  mynorm(x0./sc);            # normalement c'est la bonne valeur cf page 168
  #println("d0 = ", d0)
  d1 =  mynorm(k0./sc);              # normalement c'est la bonne valeur cf page 168
  h0 = 0.01*(d0/d1);
  if (d0 < 1.e-5) || (d1 < 1.e-5)
    h0 = 1.e-6;
  end;  
  x1 = x0 + h0*k0;
  k1 = rhs(x1,par,t0+h0);
  d2 =  mynorm((k1-k0)./sc)/h0;      # normalement c'est la bonne valeur cf page 168
  if max(d1,d2) < 1.e-15
    h1 = max(1.e-6,h0*1.e-3);
  else
    h1 = (0.01/max(d1,d2))^(1/4);  
  end
  h = min(100*h0,h1);
  return h
end
#
# Initialisation
  #Scal_Type=eltype(x0)
  p=4; # ordre
  t0=t0tf[1]; tf=t0tf[2];
  hmax=tf-t0;
  Npasmax=1000;
  n = length(x0);
  T = [t0]; 
  #T = [t0];
  #X = [Scal_Type.(x0)]
  X = [x0]
  #
  #
  # step initialisation
  h=inith(rhs,x0,par,t0,Atol,Rtol)
  nstep=0;
  t=t0
  x=x0; fin=0;

  while (nstep < Npasmax) && (fin == 0)
    nstep = nstep+1
    # Runge-Kutta method
    k1 = rhs(x,par,t);
    k2 = rhs(x+(h/3)*k1,par,t+h/3,);
    k3 = rhs(x+h*(-k1/3+k2),par,t+2*h/3);
    k4 = rhs(x+h*(k1-k2+k3),par,t+h);
    x1 = x+(h/8)*(k1+3*k2+3*k3+k4);
    xhat1 = x+(h/12)*(k1+6*k2+3*k3+2*rhs(x1,par,t));
    # err
    sc = Atol .+ max.(abs.(x),abs.(x1)) .* Rtol;
    err = mynorm((x1-xhat1)./sc);
    # calcul du pas
    if (err < 1)
      t = t+h; x = x1;
      push!(T, t)
      push!(X, x1)
      if (t > tf*(1. - eps()))
        fin = 1;
      end;  
    end;  
    h = h*min(5,max(0.2,0.9*(1/err)^(1/p)));
    if (t+h > tf)
      h = tf-t;
    end;  
  end;

    return T,X
  end


  # Fixed steps
  function myode43(rhs,x0,par,T)
    Scal_Type=eltype(x0)
    N = length(T)-1
    X = [Scal_Type.(x0)]
    x = x0
   
    for i in 1:N
      t = T[i]; h = T[i+1]-T[i]
      k1 = rhs(x,par,t);
      k2 = rhs(x+(h/3)*k1,par,t+h/3,);
      k3 = rhs(x+h*(-k1/3+k2),par,t+2*h/3);
      k4 = rhs(x+h*(k1-k2+k3),par,t+h);
      x1 = x+(h/8)*(k1+3*k2+3*k3+k4);
      x = x1;
      push!(X, x1)
    end  
    return T,X
  end

  function get_Xij(X,i,j)
    """ 
        Get from the vector of vector X the i,j composante for all times
    """
    N = length(X)
    Xij = zeros(N)
    for l in 1:N
      Xij[l] = X[l][i,j]
    end
    return Xij
  end


