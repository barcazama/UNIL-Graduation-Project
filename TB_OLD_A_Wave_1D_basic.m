%% Computing of seismic waves in 1D without boundary
clear all, close all, clc
% Variables
lenght = 1e4; % dimension of main material [m]
nx = 100; % number of mathematical boxes (augment precision and computing time) [-]
E = 1e10; % Young's modulus of material (rock) [Pa]
% E_b = 1; % High Young's modulus to kill wave in the extended space outside plot [Pa]
rho = 4500; % Density of material [kg.m⁻³]

% Formulas
B = lenght*0.1; % Extension of space to kill waves outside of the plot [m], fixed at 10% of H after and before the ploting area.
dx = (lenght+2*B)/(nx-1); % Increment on x axis [m]. We may need to add ./ for multiple dimension! [m]
dt = 1/10 * min(dx)/sqrt(E/rho); % increment of time, adapted to fit the input values [s]

% That part has to be abandonned because we cannot absorbe the wave in such
% fashion.
% E(1,1:nx-1) = E_b; % Input E_B value on the E array [Pa]
% E(1,1+nx*0.1:nx*0.9) = E_a; % Input of material's E [Pa] Should it be E(1,1+nx*0.1:(nx-1)*0.9)? Check later

X = -B:dx:lenght+B; % Creation of x array [m]

% Possible problem, the added area make it that the plot area start after 0
% because of the displacement of dx, keeping that in mind in case of futur
% problem. Workaround would be to fix the case or something more elegant as
% this will be a pain to write the position, plot like the following line:
% % plot(x(1+0.1*nx:nx-0.1*nx),E(1,1+nx*0.1:nx-0.1*nx))
S = zeros(1,nx-1); % Creation of sigma (stress) array [N.m⁻²]
V = exp(-(X-lenght/2).^2/2/(0.05*(lenght+2*B))^2); % GAUSSIAN FUNCTION BUT MUST FOUND ANOTHER WAY TO START THE WAVE


%% Plot
for i=1:600
    S = S + E.*diff(V)/dx *dt;
    V(2:nx-1) = V(2:nx-1) + 1/rho * diff(S)/dx *dt;
    plot(X,V,'-')
    axis([-B lenght+B -1 1])
    drawnow
end