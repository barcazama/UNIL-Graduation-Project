clearvars; close all; clc;

for ii = 1:2
%% Input
% Physics input
Lx = 1; % medium lenght
Ly = 1;
E = 1e1; % young's modulus
nu = 0.25; % poisson modulus
mu = 0.25; % viscosity modulus
% mu1 = 0.1; % viscosity modulus
% mu2 = 0.1; % viscosity border
rho = 1; % density
Vmax = 2; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
ny = 100;
if ii==1
    nt = 1; % number of time steps [-]
elseif ii==2
    nt = 700;
elseif ii==3
    nt = 2700;
else
    disp('Error')
end
modulo = 10; % set modulo for faster plotting
CFL = 0.1; % Courant–Friedrichs–Lewy condition making dt smaller

% mu = zeros(nx-1,ny-1);
% mu(1:round(Lx/dx),1:round(Ly/dy))=mu(1:round(Lx/dx),1:round(Ly/dy))+mu2;
% mu(round(2*Lx/dx):end,round(2*Lx/dx):end)=mu(round(2*Lx/dx):end,round(2*Lx/dx):end)+mu2;
% mu(round(Lx/dx):round(2*Lx/dx),round(Lx/dx):round(2*Lx/dx))=...
%     mu(round(Lx/dx):round(2*Lx/dx),round(Lx/dx):round(2*Lx/dx))+mu1;

%% Preprocessing
% Modeling formulas
dx = 3*Lx/(nx-1); % set dx size
dy = 3*Ly/(ny-1);
dtE = min(dx,dy)/sqrt(E/rho)*CFL;
dtV = dx^2/(mu/rho)/2*CFL;
dt = min(dtE,dtV);
x = 0:dx:3*Lx; % create position array [m]
y = 0:dy:3*Ly;
Vx = zeros(nx,ny); % displacement speed array
Vy = zeros(nx,ny);
xc = dx/2:dx:3*Lx-dx/2; % center of x axis
yc = (dy/2:dy:3*Ly-dy/2)'; % center of y axis
xc2D = repmat(xc,[nx-1,1]); % array to matrix
yc2D = repmat(yc,[1,ny-1]);
sigmaE_x = Vmax*ones (nx-1,ny-1).*exp(-((xc2D-3*Lx/2).^2+(yc2D-3*Ly/2).^2)/1/4/1e-3); %  % create intial stress Gaussian
sigmaE_y = sigmaE_x; % stress on y
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
    
    sigmaE_x = sigmaE_x +E/(1-nu*nu)*(edot_xx+nu*edot_yy)*dt;
    sigmaE_y = sigmaE_y +E/(1-nu*nu)*(edot_yy+nu*edot_xx)*dt;
    sigmaV_x = mu*(edot_xx+nu*edot_yy);
    sigmaV_y = mu*(edot_yy+nu*edot_xx);
    sigma_x = sigmaE_x+sigmaV_x;
    sigma_y = sigmaE_y+sigmaV_y;
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
        axis([Lx 2*Lx Ly 2*Ly 0 Vmax/10])
        title(['Time = ' num2str(time)])
        xlabel('[m]')
        ylabel('[m]')
        zlabel('du/dt [m.s⁻¹]')
        drawnow
    end
end
end