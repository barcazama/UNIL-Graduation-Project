clearvars; close all; clc;

%% Input
% Physics input
Lx = 1; % medium lenght
E = 1; % young's modulus
mu = 0.01; % viscosity
rho = 1; % density
Vmax = 1; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
nt = 1000000; % number of time steps [-]
CFL = 0.99; % Courant–Friedrichs–Lewy condition making dt smaller
modulo = 1000; % define modulo for fewer plot

%% Preprocessing
% Modeling formulas
dx = Lx/(nx-1); % set dx size
x = -Lx/2:dx:Lx/2; % create position array
Vx = Vmax*exp(-(x*2*pi*2).^2); % create intial displacement speed Gaussian
Vx0 = Vx; % store initial displacement speed
dt1    =   dx/(sqrt(E/rho))/2*CFL; % set dt relatif to elasticity
dt2    =   dx^2/(mu/rho)/2*CFL; % set dt relative to viscosity
dt     =   1e-5; % choose lowest dt
time = 0; % create time variable

% Physics formulas
s = E*diff(Vx)/dx*dt; % strain
s0 = s; % store intial strain
sE = s; % store intial strain for viscosity change


%% Computing and ploting wave evolution
for i=1:nt
    time = time+dt; % time counter
    edot = diff(Vx)/dx; % derivative of displacement speed
    sV = mu*edot; % calculate updated viscosity strain
    sE = sE + E*edot*dt; % calculate updated elasticity strain
    s = sV+sE; % calculate total strain according to Kelvin-Voigt model
    Vx(2:end-1) =  Vx(2:end-1) + diff(s)/dx/rho*dt; % update displacement speed
    if mod(i,modulo)==1
        plot(x,Vx)
        axis([-Lx/2 Lx/2 -1 1])
    %     plot(x(1:end-1),s)
    %     axis([-Lx/2 Lx/2 0 10])
        drawnow
    end
end