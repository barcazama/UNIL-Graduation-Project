clearvars; close all; clc;

for ii = 1:2 % 1: free boundary, 2: fixed boundary
%% Input
% Physics input
Lx = 1; % medium lenght
E = 1; % young's modulus
rho = 1; % density
Vmax = 1; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
if ii == 1
    nt = 120;
elseif ii == 2
    nt = 120;
else
    disp('Error')
end
CFL = 0.99; % Courant–Friedrichs–Lewy condition making dt smaller

%% Preprocessing
% Modeling formulas
dx = Lx/(nx-1); % set dx size
x = -Lx/2:dx:Lx/2; % create position array
Vx = Vmax*exp(-(x*2*pi*2).^2); % create intial displacement speed Gaussian
Vx0 = Vx; % store initial displacement speed
s = zeros(1,nx-1); % create stress array
dt    =   dx/(sqrt(E/rho))/2*CFL; % set dt relatif to elasticity
time = 0; % create time variable

% Physics formulas



%% Computing and ploting wave evolution
figure(ii)
for i=1:nt
    time = time+dt; % time counter
    if ii == 1 % free boundary
        s(1)=0;
        s(end)=0;
    elseif ii == 2 % fixed boundary
        s(1)=s(2);
        s(end) = s(end-1);
    end
    s = s + E*diff(Vx)/dx*dt;
    Vx(2:nx-1) = Vx(2:nx-1) +diff(s)/dx/rho*dt;
    plot(x,Vx)
    axis([-Lx/2 Lx/2 -Vmax Vmax])
    xlabel('Position [m]')
    ylabel('Displacement speed [m.s⁻¹]')
    if ii == 1 % free boundary
        title('Model with free boundary')
        s(1)=0;
        s(end)=0;
    elseif ii == 2 % fixed boundary
        title('Model with fixed boundary')
        s(1)=s(2);
        s(end) = s(end-1);
    end
    grid on; box on;
    drawnow
end

end