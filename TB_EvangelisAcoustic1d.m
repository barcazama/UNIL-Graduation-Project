clear,clc,figure(1),clf
%--------------
%Physics
Maxwell =  0; % 1 = maxwell equation, else kelvin
Lx     =   1; % medium lenght on x axis
E      =   1; % youngs modulus
eta    =   1; % viscosity
rho    =   1; % density
Vlam   =   1;
V0     =   0; 
%Numerics
nx     =   100; % number of numerical cases
CFL    =   0.99; % factor making dt 1% smaller
nout   =   100; % modulo factor
nt     =   500*nout; % number of plot
%Preprocess
dx     =   Lx/(nx-1);
x1     =   -Lx/2:dx:Lx/2;
x1e    =   [x1,x1(end)+dx/2];
Vx     =   V0 + Vlam*exp(-(x1e*2*pi*2).^2);
dt1    =   dx/(sqrt(E/rho))/2*CFL;
dt2    =   dx^2/(eta/rho)/2*CFL;
dt     =   min(dt1,dt2);
s      =   dt*E*diff(Vx)/dx;
%s      =   0 + 1*exp(-(x1*2*pi*2).^2); Vx = Vx*0;
s0     =   s;
sE     =   s;
Vx0    =   Vx;
time   =   0;
strain =   0;
for it = 1:nt
    time        =  time+dt;
    edot        =  diff(Vx)/dx;
    if Maxwell == 1
        ksi         =  eta/dt/E;
        s           =  ksi*s/(1+ksi) + eta*edot/(1+ksi);
    else %Kelvin
        sV          =  eta*edot;
        sE          =  sE + dt*E*edot;
        s           =  sV + sE;
    end    
    strain      =  strain + edot*dt;
    Vx(1)       =  Vx(2);
    Vx(end)     =  Vx(end-1);
    Vx(2:end-1) =  Vx(2:end-1) + dt*diff(s)/dx/rho;
          
    if mod(it,nout)==1
        subplot(121)
        plot(x1,s,x1,s0)
        axis([-Lx/2 Lx/2 -0.5 1])
        subplot(122)
        %plot(x1e,Vx,x1e,Vx0)
        %axis([-Lx/2 Lx/2 -6 6])
        plot(strain,s,'o')
        axis([0 1 0 1]),axis square
        drawnow
    end    
end

