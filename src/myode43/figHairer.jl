#
# ~gergaud/ENS/Control/ODE/pasvar/figHairer.m
#
# Auteurs:  Joseph GERGAUD
# Date:     fev. 2007
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2, rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# pour visualiser les graphiques du livre 
# ref: Hairer page 170 T1
# Nombre de pas acceptes et rejetes pas variables
include("myode43_for_graph.jl")
using Plots

function bruss(x, par, t)
    λ = par[1]
    x₁ = x[1]
    x₂ = x[2]
    return [1+x₁^2*x₂-(λ+1)*x₁, λ*x₁-x₁^2*x₂]
end

t0=0;
tf=20;
# 
Atol=1.e-4;
Rtol=1.e-4;
x0=[1.5, 3];
par = [3.0]
T, X, haccept, Trej, hrej, Ttot, err_loc_est, err_loc_exact, err_abs = myode43(
    bruss, x0, par, [t0, tf], Rtol, Atol
)

N = length(T)
println("N = ", N)
println(size(X))

X1 = [X[i][1] for i in 1:N]
p1 = plot(T, X1; labels=false)

scatter!(p1, T, X1; c=:blue, labels=false)
X2 = [X[i][2] for i in 1:N]
plot!(p1, T, X2; c=:green, labels=false)
scatter!(p1, T, X2; c=:green, labels=false, ylabel="x_1, x_2")
#

npasaccept = length(haccept);
p2 = plot(T[1:npasaccept], haccept; yscale=:log10, labels=false);
scatter!(p2, T[1:npasaccept], haccept; yscale=:log10, c=:blue, labels="accepted steps");
scatter!(p2, Trej, hrej; yscale=:log10, c=:red, labels="reject steps");

p3 = plot(Ttot[2:end], err_loc_est; yscale=:log10, c=:blue, labels="local estimated errors")
scatter!(
    p3,
    Ttot[2:end],
    err_loc_est;
    yscale=:log10,
    c=:blue,
    labels=false,
    markersize=3,
    markershape=:x,
)
plot!(p3, [t0, tf], [Atol, Atol]; c=:black, labels=false)

plot!(
    p3, Ttot[2:end], err_loc_exact; c=:green, yscale=:log10, labels="local estimated errors"
)
scatter!(
    p3,
    Ttot[2:end],
    err_loc_exact;
    c=:green,
    yscale=:log10,
    labels=false,
    markersize=3,
    markershape=:x,
)
plot!(p3, Ttot[2:end], err_abs; c=:black, yscale=:log10, labels="global exact errors")
scatter!(
    p3,
    Ttot[2:end],
    err_abs;
    c=:black,
    yscale=:log10,
    labels=false,
    markershape=:x,
    markersize=3,
)

plot(p1, p2, p3; layout=(3, 1))

#=
#
# Equation avec discintinuit� du second membre
# Figure page 197 avec mon programme pas variable
t0=0;tf=1;
# 
Atol=0.7e-3; Rtol=0.7e-3;
y0=0.3;
[T,Y,haccept,Trej,hrej,Ttot,err_loc_est,err_loc_exact,err_abs]=myode43(@phi_disc,[t0 tf],y0,Rtol,Atol);
figure;
j = 1;
for i=1:length(Ttot),
  if T(j)==Ttot(i),
    I(i)=1;
    j=j+1;
  else
    I(i)=0;   
  end;   
end;
h1=plot(Ttot(find(I==1)),find(I==1)-1,'o')
hold on
h2=plot(Ttot(find(I==0)),find(I==0)-1,'xr')
plot(Ttot,0:length(Ttot)-1)
xlabel('t');
ylabel('Nombre de pas');
legend([h1 h2],'pas accept�s','pas rejet�s')
print('fig_disc1','-depsc')

=#
