clear all, close all, clc

% nx = 11;
% width = 1;
% X = 0:width/(nx-1):1;
% dx = X(2)-X(1);
% kappa = 1;
% T_old = sin(X/width*pi);
% T_new = T_old;
% dt = 1e-3;
% no_dt = 100;
% Ind = 2:nx-1;
% 
% for t=1:no_dt
%    T_new(Ind) = T_old(Ind) + dt*kappa/dx^2 * (T_old(Ind+1)-2*T_old(Ind)+T_old(Ind-1));
%    T_old(Ind) = T_new(Ind);
%    plot(X,T_new,'-o')
%    axis([o width o 1])
%    xlabel('X','fontsize',5)
%    ylabel('T','fontsize',5)
%    drawnow
% end

%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1;
E = 1e1;
nu = 0.25;
H = [1,1];
nx = [71,71];


dx = H./(nx-1);
dt = 1/10*min(dx)/sqrt(E/rho);
x = [0:dx(1):H(1)];
y = [0:dx(2):H(2)];
x2D = repmat(x,[1,nx(2)]);
y2D = repmat(y,[nx(1),1]);

Vx = zeros(nx(1),nx(2));
Vy = zeros(nx(1),nx(2));
xc = [dx(1)/2:dx(1):H(1)-dx(1)/2]';
yc = [dx(2)/2:dx(2):H(2)-dx(2)/2];
xc2D = repmat(xc,[1,nx(2)-1]);
yc2D = repmat(yc,[nx(1)-1,1]);
Sxx = ones (nx(1)-1,nx(2)-1).*3e-3/1e-3.*exp(-((xc2D-H(1)/2).^2+(yc2D-H(2)/2).^2)/1/4/1e-3);
Txy = zeros(nx(1)-1,nx(2)-1);
Syy = Sxx;

%%

for it = 1:710
    dVx_dx = (diff(Vx(:,1:nx(2)-1),1,1)+diff(Vx(:,2:nx(2)),1,1))/dx(1)/2;
    dVy_dx = (diff(Vy(:,1:nx(2)-1),1,1)+diff(Vy(:,2:nx(2)),1,1))/dx(1)/2;
    dVx_dy = (diff(Vx(1:nx(1)-1,:),1,2)+diff(Vx(2:nx(1),:),1,2))/dx(2)/2;
    dVy_dy = (diff(Vy(1:nx(1)-1,:),1,2)+diff(Vy(2:nx(1),:),1,2))/dx(2)/2;
    
    Sxx = Sxx +dt*E/(1-nu*nu)*(dVx_dx+nu*dVy_dy);
    Syy = Syy +dt*E/(1-nu*nu)*(dVy_dy+nu*dVx_dx);
    Txy = Txy +dt*E/(1+nu)/2*(dVx_dy+dVy_dx);
    
    D_Sxx_x = (diff(Sxx(:,1:nx(2)-2),1,1)+diff(Sxx(:,2:nx(2)-1),1,1))/dx(1)/2;
    D_Txy_x = (diff(Txy(:,1:nx(2)-2),1,1)+diff(Txy(:,2:nx(2)-1),1,1))/dx(1)/2;
    D_Txy_y = (diff(Txy(1:nx(1)-2,:),1,2)+diff(Txy(2:nx(1)-1,:),1,2))/dx(2)/2;
    D_Syy_y = (diff(Syy(1:nx(1)-2,:),1,2)+diff(Syy(2:nx(1)-1,:),1,2))/dx(2)/2;
    
    Vx(2:nx(1)-1,2:nx(2)-1) = Vx(2:nx(1)-1,2:nx(2)-1)+dt*(D_Sxx_x+D_Txy_y)/rho;
    Vy(2:nx(1)-1,2:nx(2)-1) = Vy(2:nx(1)-1,2:nx(2)-1)+dt*(D_Txy_x+D_Syy_y)/rho;
    
    if mod(it,10)==1
        surf(x,y,(Vx'.^2+Vy'.^2).^(1/2)); axis equal; axis tight; shading interp; light; lighting phong;
        view(30,40),caxis([0 0.075]);axis([0 1 0 1 0 0.2])
        drawnow
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





