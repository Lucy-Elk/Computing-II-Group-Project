% What this script does:
% Uses u and w from 'wind_turbine_results.m'
% 'Seeds' starting positions of particles
% Computes how each particle moves over time, using particle_trajectory.m
% Plots the trajectories, x against z

% Inputs:
% u,w - velocities from (wind_turbine_results.m)

% Outputs:
% Plots of trajectories

clear; clc;


% load results from tasks 4 and 5
load('wind_turbine_results.m');

%Set up:
Lx=x(end)-x(1);
Lz= z(end)-z(1);

% ChatGPT Fix
x_min=x(1);
z_min=z(1);

x_turb= x_0; % turbine x position
z_hub= z_0; % turbine hub height

% 'Seed' the starting points
x_seeds = linspace(x_min + 0.1*Lx, x_min + 0.4*Lx, 5);
z_seeds = [50,100,150,200,250,300,350,400,450,500];

[x_seeds_grid,z_seeds_grid]= meshgrid(x_seeds,z_seeds);
X_seeds= x_seeds_grid(:);
Z_seeds = z_seeds_grid(:);
numP= numel(X_seeds);

% Initial setup
tmax=(Lx/U0)*1.5;
dt= tmax/1000;
nt= round(tmax/dt)+1;
xp_all=zeros(numP,nt);
zp_all=zeros(numP,nt);

% Loop over particles
for p=1:numP
    x0=X_seeds(p);
    z0=Z_seeds(p);
    
    % Use particle_trajectory.m to find trajectories
    [xp,zp,t]=particle_trajectory(x0,z0,x,z,u,w,U0,dt,tmax);
    
    xp_all(p,:) = xp; 
    zp_all(p,:) = zp; 
end

% Plot trajectories

figure;
hold on;

% Plot particle trajectories

colors = lines(length(z_seeds));

% for p = 1:numP
% 
%     plot(xp_all(p,:), zp_all(p,:), 'LineWidth', 1.2);
% end

for p = 1:numP
    % Copy particle trajectory
    xp_plot = xp_all(p,:);
    zp_plot = zp_all(p,:);
    
    % unwrap for visual ease
    for n = 2:nt
        if xp_plot(n) < xp_plot(n-1) - 0.5*Lx
            xp_plot(n:end) = xp_plot(n:end) + Lx;
        elseif xp_plot(n) > xp_plot(n-1) + 0.5*Lx
            xp_plot(n:end) = xp_plot(n:end) - Lx;
        end
    end
    % get which z level the particle came from
    z_level = find(abs(Z_seeds(p) - z_seeds) < 1e-6);
    plot(xp_plot, zp_plot, 'LineWidth', 1, 'Color',colors(z_level,:));
end

% Mark where turbine is
plot(x_turb,z_hub,'rx','Markersize',10,'Linewidth',2);
text(x_turb, z_hub, '  Turbine', ...
     'Color', 'r', 'FontSize', 20, 'FontWeight', 'bold', ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');


xlabel('X Position (m)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Z Position (m)', 'Interpreter', 'latex', 'FontSize', 24);
title('Particle Trajectories', 'Interpreter', 'latex', 'FontSize', 26);
xlim([x(1),x(end)]);
ylim([z(1),525]);
set(gca, 'FontSize', 24, 'TickLabelInterpreter', 'latex');
grid on;
hold off;



