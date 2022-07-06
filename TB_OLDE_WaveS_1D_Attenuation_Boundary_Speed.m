%% Computing of seismic waves in 1D without boundary
clear all, close all, clc
% Variables
lenght = 5e4; % dimension of main material [m]
nx = 1e4; % number of mathematical boxes (augment precision and computing time) [-]
E = 1e10; % Young's modulus of material (rock) [Pa]
rho = 4500; % Density of material [kg.m⁻³]
Qp_0 = 1; % Anelastic attenuation factor of layer 0 and P wave[-]
Qs_0 = 1; % Anelastic attenuation factor of layer 0 and S wave[-]
Vpx = 6000; % Horizontal speed of P wave[m.s⁻¹]
Vsx = 2000;

% Formulas
B = lenght*0.4; % Extension of space to kill waves outside of the plot [m], fixed at 40% of H after and before the ploting area.
dx = (lenght+2*B)/(nx-1); % Increment on x axis [m]. We may need to add ./ for multiple dimension! [m]
dt = 1/2.1 * min(dx)/sqrt(E/rho); % increment of time, adapted to fit the input values, especially dx [s]
% dt = (lenght+2*B)/Vpx/nx; % number of second by mathematical boxes to match wave speed [s]

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
Vu = exp(-(X-lenght/2).^2/2/(0.05*(lenght+2*B))^2); % Initialise a starting velocity based on a Gaussian. Vu reprensent the displacement perdicular of the wave direction [m.s⁻¹]

%% Creation of Q array to adapt Q in boundary
Qp = zeros(1,nx); % Create an empty array for anelastic attenuation factor [-]
for i=1:nx % Loop to reduce Q in the extended area
    if X(i)<0 % Target area before plot
        Qp(i)=Qp_0-Qp_0*erf(-X(i)/B); % Slowly reduce Q in that area using the error function 
    elseif X(i)>lenght % Target area after plot
        Qp(i)=Qp_0-Qp_0*erf((X(i)-lenght)/B); % Slowly reduce Q in that area using error function
    else
        Qp(i) = Qp_0; % Assign Q_0 value in the initial layer
    end
end

%% Calcul of Period
S_1 = S; % Create a new S_1 array similar to S for manipulation without affecting the plot [N.m⁻²]
Vu_1 = Vu; % Same but for V [m.s⁻¹]
c = round(nx/2); % Stock the middle position ofQ = 1; % Anelastic attenuation factor [-] the array [-]
Tp = 0; % Create empty period value [s]

while Vu_1(c)>0
    Vu_0 = Vu_1; % Stock previous value of V_1 [m.s⁻¹]
    S_1 = S_1 + E*diff(Vu_1)/dx*dt; % Calculate new itération sigma [N.m⁻²]
    Vu_1(2:nx-1) = Vu_1(2:nx-1) + 1/rho * diff(S_1)/dx *dt; % Calculate new itération displacement speed [m.s⁻¹]
    Tp = Tp+dt; % Calculate half period of P wave [s]
end
Tp = Tp*2; % Change half period to period of P wave [s]
    
%% Plot
Vu_1 = Vu; % Reset V_1 to V values [m.s⁻¹]
i = 0; % Creation of itération variable [-]
r = 1e9; % Creation of break variable to initialize while loop [-]
time = 0; % Creation of time variable [s]
while i<1000000
    time = time + dt; % Counting time [s]
    i = i + 1; % Update number of itération [-]
    Vu_0 = Vu_1; % Stock previous itération V_1 value (used to calcul the break function) [m.s⁻¹]
    S = S + E.*diff(Vu_1)/dx *dt; % Calculate new itération sigma [N.m⁻²]
    Vu_1(2:nx-1) = Vu_1(2:nx-1) + 1/rho * diff(S)/dx *dt ;
    
%     Vu_1(2:nx-1) = sign(Vu_1(2:nx-1)).*(Vu_1(2:nx-1).^2-2*pi.*Vu_1(2:nx-1).^2./Qp(2:nx-1).*dt/Tp).^(1/2); % Calculate new itération displacement speed [m.s⁻¹]
    Vu_1(2:nx-1)=Vu_1(2:nx-1).*(1-(1-(1-2*pi./Qp(2:nx-1)).^(1/2))*dt/Tp);
    r = sum(abs(Vu_1-Vu_0)); % Calculating new trigger break condition, rest between V_1 and V_0

% Formula
    if mod(i,600)==1
        plot(X,Vu_1,'-')
        title([num2str(time),' s'])
        axis([-B lenght+B -1 1])
        drawnow
    end
end