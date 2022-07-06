clearvars; close all; clc;

for ii = 1:3
%% Input
% Physics input
Lx = 1; % medium lenght
E = 1; % young's modulus
rho = 1; % density
Vmax = 1; % max value of the initial Gaussian function

% Numericals input
nx = 100; % number of grid points [-]
if ii == 1
    nt = 60;
elseif ii == 2
    nt = 150;
elseif ii ==3
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
    if ii == 3
        s(1)=0;
        s(end)=0;
    end
    s = s + E*diff(Vx)/dx*dt;
    Vx(2:nx-1) = Vx(2:nx-1) +diff(s)/dx/rho*dt;
    plot(x,Vx)
    axis([-Lx/2 Lx/2 -Vmax Vmax])
    title(['Time = ' num2str(time)])
    xlabel('Position [m]')
    ylabel('Displacement speed [m.s⁻¹]')
    grid on; box on;
    drawnow
end

end