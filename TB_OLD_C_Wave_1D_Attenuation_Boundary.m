%% Computing of seismic waves in 1D without boundary
clear all, close all, clc
% Variables
lenght = 5e4; % dimension of main material [m]
nx = 1e4; % number of mathematical boxes (augment precision and computing time) [-]
E = 1e10; % Young's modulus of material (rock) [Pa]
rho = 4500; % Density of material [kg.m⁻³]
Q_0 = 1; % Anelastic attenuation factor [-]

% Formulas
B = lenght*0.2; % Extension of space to kill waves outside of the plot [m], fixed at 20% of H after and before the ploting area.
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

%% Creation of Q array to adapt Q in boundary
Q = zeros(1,nx); % Create an empty array for anelastic attenuation factor [-]
for i=1:nx
    if X(i)<0 % Target area before plot
        Q(i)=Q_0-Q_0*erf(-X(i)/B); % Slowly reduce Q in that area using the error function 
    elseif X(i)>lenght % Target area after plot
        Q(i)=Q_0-Q_0*erf((X(i)-lenght)/B); % Slowly reduce Q in that area using error function
    else
        Q(i) = Q_0; % Assign Q_0 value in the initial layer
    end
end

%% Calcul of Period
S_1 = S; % Create a new S_1 array similar to S for manipulation without affecting the plot [N.m⁻²]
V_1 = V; % Same but for V [m.s⁻¹]
c = round(nx/2); % Stock the middle position ofQ = 1; % Anelastic attenuation factor [-] the array [-]
T = 0; % Create empty period value [s]

while V_1(c)>0
    V_0 = V_1; % Stock previous value of V_1 [m.s⁻¹]
    S_1 = S_1 + E*diff(V_1)/dx*dt; % Calculate new itération sigma [N.m⁻²]
    V_1(2:nx-1) = V_1(2:nx-1) + 1/rho * diff(S_1)/dx *dt; % Calculate new itération displacement speed [m.s⁻¹]
    T = T+dt; % Calculate half period [s]
end
T = T*2; % Change half period to period [s]
    
%% Plot
V_1 = V; % Reset V_1 to V values [m.s⁻¹]
i = 0; % Creation of itération variable [-]
r = 1e9; % Creation of break variable to initialize while loop [-]

while r>3e-3
    i = i + 1; % Update number of itération [-]
    V_0 = V_1; % Stock previous itération V_1 value [m.s⁻¹]
    S = S + E.*diff(V_1)/dx *dt; % Calculate new itération sigma [N.m⁻²]
    V_1(2:nx-1) = V_1(2:nx-1) + 1/rho * diff(S)/dx *dt - V_1(2:nx-1)./Q(2:nx-1)*2*pi/T*dt; % Calculate new itération displacement speed [m.s⁻¹]
    r = sum(abs(V_1-V_0)); % Calculating new trigger break condition, rest
    
    
    if mod(i,200)==1
        plot(X,V_1,'-')
        axis([-B lenght+B -1 1])
        drawnow
    end
end