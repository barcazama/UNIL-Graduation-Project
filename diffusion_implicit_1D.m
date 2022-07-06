clear, clc, close all 
%% parametres
Nx = 40;
Nt = 1000;
dx = 1/(Nx-1);
c = 5;
l=1;
cfl = (2*c)/(2*c+1);
dt = cfl*(dx^2)/(2*c);
alpha = c*dt/dx^2;
beta = (1+2*alpha);
x = 0:dx:l;
%% construction de la matrice A
A = full(spdiags(ones(Nx,1)*[-alpha,beta,-alpha],[-1,0,1],Nx,Nx));
%% construction de la matrice b 
b = ones(Nx,1);
b(1:Nx/4,1)=0;         % initial conditions
b(3*Nx/4:Nx,1)=0;      % initial conditions
%% %% %% %%t
time = 0;
for t = 0:dt:Nt
   
    T = linsolve(A,b);
    
    T(1)     =  0;     % boundary conditions
    T(Nx)    =   1;    % boundary conditions
    b = T;
    time        =   time+dt;
    
    figure(1), clf     % Solution plot
    plot(x,T);
    xlabel('x [m]')
    ylabel('Temperature [^oC]')
    ylim([0 3])
    drawnow
end