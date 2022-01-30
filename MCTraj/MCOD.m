function [] = MCOD(nElec)


    % clear all
    
    %clearvars
    %clearvars -GLOBAL
    %close all
    
    % set(0,'DefaultFigureWindowStyle','docked')
    global C

    addpath ../geom2d/geom2d

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      %metres (32.1740 ft) per s²

    accel = 1;
    V0 = 0;
    dt = 15e-15;
    %type = 0;
    doPlot = 1;
    MarkerSize = 12;
    TStop = 5e-12;
    t = 0;
    Size = 2e-25;
    Limits = [-0.1*Size 10*Size -Size +Size];
    LimitsTime = [0 5e-12 -1e-15 10e-13];
    Vavg = 0;
    count = 1;
    
    

    x(1, :) = zeros(1, nElec);
    Vx(1:nElec) = V0;
    
    while t < TStop

            R = rand(1, nElec); % Make random array

            dvx = accel * dt; % get acceleration
            %dvx = randn(1, nElec) * paras(1)*dt;
            Vx = Vx + dvx; % get velocity
            for i = 1:nElec % check for scatter
               if R(i) <= 0.05
                   Vx(i) = Vx(i) * -0.25; % get scatter velocity
               else
               end
            end
            dx = Vx * dt; % get displacement
            x(1,:) = x(1,:) + dx; % check for position change
            Vavg = sum(Vx)/nElec; % get current average Velocity
            VavgPlot(count,:) = [Vavg]; % store array of times and average velocities
            timePlot(count,:) = [t];
            count = count + 1;
            t  = t + dt; % iterate time 
            
            
        if doPlot
        
        subplot(2,1,1); %plot the animation
            plot(x, zeros(size(x)), 'bo', 'markers',...
            MarkerSize,'MarkerFaceColor', 'b');
            hold on
            subplot(2,1,1), plot(x, zeros(size(x)), 'bo', 'markers',...
            MarkerSize,'MarkerFaceColor', 'b');
            quiver(x,zeros(size(x)),Vx ,zeros(size(x)));
            hold off
            axis(Limits);
            xlabel('x');
            ylabel('y');
            grid on
         subplot(2,1,2); %plot Avg Velocity
            plot(t,Vavg,'v','linewidth', 2);
            hold on
            subplot(2,1,2), plot(timePlot,VavgPlot,'linewidth', 2);
            hold off
            axis(LimitsTime);
            xlabel('Time');
            ylabel('Velocity');
            grid on
            
            pause(0.0001)
            
        end
    end
end

