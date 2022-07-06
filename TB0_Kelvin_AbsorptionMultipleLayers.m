clearvars; close all; clc;

for ii= 1:3
%% Input
% Physics input
Lx1 = 1; % first medium lenght
Lx2 = 1; % second medium length
Lx3 = 1; % outside boundary medium lenght
E1 = 0.01; % young's modulus
E2 = 0.02;
E3 = 1;
mu1 = 0.01; % viscosity
mu2 = 0.02; 
mu3 = 0.03; % out of boundary viscosity value
rho = 1; % density
Vmax = 1; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
if ii == 3
    nt = 1000000; % number of time steps [-]
elseif ii == 2
    nt = 500000;
elseif ii == 1
    nt = 1;
else
    disp('Error')
end
CFL = 0.99; % Courant–Friedrichs–Lewy condition making dt smaller
modulo = 5000; % define modulo for fewer plot

%% Preprocessing
% Modeling formulas
Lx = Lx1+Lx2*2+Lx3*2; % total lenght of all the medium
dx = Lx/(nx-1); % set dx size

% mu1 = zeros(1,Lx1/dx)+mu1; % creating array with relative viscosity for each medium
% mu2 = zeros(1,Lx2/dx/2)+mu2;
% mu3 = zeros(1,Lx3/dx/2)+mu3;
% mu = [mu3 mu2 mu1 mu3];
mu = zeros(1,nx-1);
mu(1+round((Lx2+Lx3)/dx):end-round((Lx2+Lx3)/dx)) = ...
    mu(1+round((Lx2+Lx3)/dx):end-round((Lx2+Lx3)/dx))+mu1;
mu(round(Lx3/dx):round(Lx3/dx)+round(Lx2/dx)) = ...
    mu(round(Lx3/dx):round(Lx3/dx)+round(Lx2/dx))...
    +mu2;
mu(end-round(Lx3/dx)-round(Lx2/dx):end-round(Lx3/dx)) = ...
    mu(end-round(Lx3/dx)-round(Lx2/dx):end-round(Lx3/dx))...
    +mu2;
mu(1:round(Lx3/dx)) = mu(1:round(Lx3/dx))+mu3;
mu(end-round(Lx3/dx):end) = mu(end-round(Lx3/dx):end)+mu3;

% E1 = zeros(1,Lx1/dx)+mu1; % creating array with relative young for each medium
% E2 = zeros(1,Lx2/dx)+mu2;
% E3 = zeros(1,Lx3/dx)+mu3;
% E = [E1 E2 E3];
E = E1;

x = -Lx/2:dx:Lx/2; % create position array
Vx = Vmax*exp(-(x*2*pi*2).^2); % create intial displacement speed Gaussian
Vx0 = Vx; % store initial displacement speed
dt1    =   dx/(sqrt(max(E)/rho))/2*CFL; % set dt relatif to elasticity
dt2    =   dx^2/(max(mu)/rho)/2*CFL; % set dt relative to viscosity
dt     =   1e-5; % choose lowest dt
time = 0; % create time variable


% Physics formulas
s = E.*diff(Vx)/dx*dt; % strain
s0 = s; % store intial strain
sE = s; % store intial strain for viscosity change


%% Computing and ploting wave evolution
figure(ii)
for i=1:nt
    time = time+dt; % time counter
    edot = diff(Vx)/dx; % derivative of displacement speed
    sV = mu.*edot; % calculate updated viscosity strain
    sE = sE + E.*edot*dt; % calculate updated elasticity strain
    s = sV+sE; % calculate total strain according to Kelvin-Voigt model
    Vx(2:end-1) =  Vx(2:end-1) + diff(s)/dx/rho*dt; % update displacement speed

    if mod(i,modulo)==1
%         subplot(211)
        plot(x,Vx,[Lx2/2 Lx2/2],[-Vmax Vmax],[-Lx2/2 -Lx2/2],[-Vmax Vmax])
        axis([-(Lx1+Lx2)/2 (Lx1+Lx2)/2 -Vmax Vmax])
        title(['Time = ' num2str(time)])
        xlabel('Position [m]')
        ylabel('Displacement speed [m.s⁻¹]')
        grid on; box on;
%         subplot(212)
%         plot(x(1:end-1),s)
%         axis([-(Lx1+Lx2)/2 (Lx1+Lx2)/2 -Vmax Vmax])
        drawnow
    end
end
end