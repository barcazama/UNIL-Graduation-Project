clearvars; close all; clc;

for ii = 1:3
%% Input
% Physics input
Lx = 1; % medium lenght
Ly = 1;
E = 1e1; % young's modulus
nu = 0.25; % poisson modulus
rho = 1; % density
Vmax = 2; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
ny = 100;
if ii==1
    nt = 1; % number of time steps [-]
elseif ii==2
    nt = 200;
elseif ii==3
    nt = 500;
end
modulo = 10; % set modulo for faster plotting
CFL = 0.1; % Courant–Friedrichs–Lewy condition making dt smaller

%% Preprocessing
% Modeling formulas
dx = Lx/(nx-1); % set dx size
dy = Ly/(ny-1);
dt = min(dx,dy)/sqrt(E/rho)*CFL; % 
x = 0:dx:Lx; % create position array [m]
y = 0:dy:Ly;
Vx = zeros(nx,ny); % displacement speed array
Vy = zeros(nx,ny);
xc = dx/2:dx:Lx-dx/2; % center of x axis
yc = (dy/2:dy:Ly-dy/2)'; % center of y axis
xc2D = repmat(xc,[nx-1,1]); % array to matrix
yc2D = repmat(yc,[1,ny-1]);
sigma_x = Vmax*ones (nx-1,ny-1).*exp(-((xc2D-Lx/2).^2+(yc2D-Ly/2).^2)/1/4/1e-3); %  % create intial stress Gaussian
sigma_y = sigma_x; % stress on y
tau_z = zeros(nx-1,ny-1); % shear stress
time = 0;

%% Computing and ploting wave evolution
figure(ii)
for i = 1:nt
    time = time+dt;
    edot_xx = (diff(Vx(:,1:ny-1),1,1)+diff(Vx(:,2:ny),1,1))/dx/2; % derivative dVx/dx derivative
    edot_yx = (diff(Vy(:,1:ny-1),1,1)+diff(Vy(:,2:ny),1,1))/dx/2;
    edot_xy = (diff(Vx(1:nx-1,:),1,2)+diff(Vx(2:nx,:),1,2))/dy/2;
    edot_yy = (diff(Vy(1:nx-1,:),1,2)+diff(Vy(2:nx,:),1,2))/dy/2;
    
    sigma_x = sigma_x +E/(1-nu*nu)*(edot_xx+nu*edot_yy)*dt;
    sigma_y = sigma_y +E/(1-nu*nu)*(edot_yy+nu*edot_xx)*dt;
    tau_z = tau_z +E/(1+nu)/2*(edot_xy+edot_yx)*dt;
    
    sdot_xx = (diff(sigma_x(:,1:ny-2),1,1)+diff(sigma_x(:,2:ny-1),1,1))/dx/2; % derivative of dsigma/dx
    tdot_x = (diff(tau_z(:,1:ny-2),1,1)+diff(tau_z(:,2:ny-1),1,1))/dx/2;
    tdot_y = (diff(tau_z(1:nx-2,:),1,2)+diff(tau_z(2:nx-1,:),1,2))/dy/2;
    sdot_yy = (diff(sigma_y(1:nx-2,:),1,2)+diff(sigma_y(2:nx-1,:),1,2))/dy/2;
    
    Vx(2:nx-1,2:ny-1) = Vx(2:nx-1,2:ny-1)+(sdot_xx+tdot_y)/rho*dt;
    Vy(2:nx-1,2:ny-1) = Vy(2:nx-1,2:ny-1)+(tdot_x+sdot_yy)/rho*dt;
    
    if mod(i,modulo)==1
        mesh(x,y,(Vx'.^2+Vy'.^2).^(1/2))
        axis equal; axis tight;
        axis([x(1,1) x(1,end) y(1,1) y(1,end) 0 Vmax/10])
        title(['Time = ' num2str(time)])
        xlabel('[m]')
        ylabel('[m]')
        zlabel('du/dt [m.s⁻¹]')
        drawnow
    end
end
end