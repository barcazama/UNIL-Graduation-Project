
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       FDTD.m
%
%    Author: M.V. Barnhart
%      Date: July 2016
% 
%   Description: Finite-difference time-domain method for evaluating
%                wave propagation in 2D.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nx = 500; % x-data points
    Ny = 500; % y-data points
    Nt = 500; % Total number of time steps
    dt = 0.5; % Time step-size

    t = 0:dt:Nt; % Time array

    u = zeros(Nx, Ny); % Initial displacement field u(x, y)
    u_p = zeros(Nx, Ny);
    u_n = zeros(Nx, Ny);

    % Gaussian source variables
    t0 = 10;
    tp = 5;

    for k = 1:length(t)

        fprintf('Calculating displacement fields: %.2f %%',k/length(t)*100)

        for i = 2:Nx-1

            for j = 2:Ny-1    

                u_n(i,j) = 2*u_n(i,j) - u_p(i,j) + ...
                          0.5*( u(i+1,j) + u(i-1,j) - ...
                          4*u(i,j) + u(i,j+1) + u(i,j-1) );

            end

        end

        u_p = u;   % Set present displacement state to previous state
        u = u_n;   % Set next displacement state to current state

        % Gaussian source at corner of grid
        %u(Nx-1,Ny-1)= exp(-((k-t0)/tp)^2);  

        % Gaussian source at center of grid
        u(Nx/2, Ny/2) = exp(-((k - t0)/tp)^2);
        
        u(Nx/4, Ny/4) = exp(-((k - t0)/tp)^2);
        u(3*Nx/4, 3*Ny/4) = exp(-((k - t0)/tp)^2);
        u(3*Nx/4, Ny/4) = -exp(-((k - t0)/tp)^2);
        u(Nx/4, 3*Ny/4) = -exp(-((k - t0)/tp)^2);
        % Initial displacement on boundaries
        %u(1:Nx,1)=0.1*real(exp(-((k-t0)/tp)^2)); % Left boundary
        %u(1:Nx,Ny)=0.1*real(exp(-((k-t0)/tp)^2)); % Right boundary

        % Plot 2D displacement field and getframe()
        surf(u);
        shading interp;colormap('jet');
        set(gcf, 'Position', get(0, 'Screensize'));
        axis([0 Nx 0 Ny -1 1]);
        axis(gca,'off');
        A(k) = getframe(gcf);
        clc;

    end

    % Generate video from frames
    disp('Generating video...');
    v = VideoWriter('newfile12','Uncompressed AVI');
    open(v);
    writeVideo(v,A);
    close(v);
    disp('Video creation complete!');
    