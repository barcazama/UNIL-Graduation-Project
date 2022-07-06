clear all, close all, clc
nx = 101;
width = 1;
X = 0:width/(nx-1):1;
dx = X(2)-X(1);
E = 1;
rho = 100;
xc = 0.05;
V = exp(-(X-width/2).^2/2/xc^2);
S = zeros(1,nx-1);
% dt = 3.3e-3;
dt = 1/10 * min(dx)/sqrt(E/rho);
no_dt = 600;
Ind = [2:nx-1];
%%
for t=1:no_dt
    S = S + E * diff(V)/dx * dt;
    V(Ind) = V(Ind) + 1/rho * diff(S)/dx *dt;
    plot(X,V,'-')
    axis([0 width -1 1])
    drawnow
end
xlabel('X','fontsize',5)
ylabel('T','fontsize',5)