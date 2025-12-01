% What this script does:
% Uses u and w from 'wind_turbine_results.m'
% 'Seeds' starting positions of particles
% Computes how each particle moves over time, using particle_trajectory.m
% Plots the trajectories, x against z

% Inputs:
% u,w - velocities from (wind_turbine_results.m)

% Outputs:
% Plots of trajectories

clear; clc; close all;


% load results from tasks 4 and 5
load('wind_turbine_results.mat');

%Set up:
Lx=x(end)-x(1);
Lz= z(end)-z(1);

% ChatGPT Fix
x_min=x(1);
z_min=z(1);

x_turb= x_turbine; % turbine x position
z_hub= z_turbine; % turbine hub height

% 'Seed' the starting points
x_seeds = linspace(x_min + 0.1*Lx, x_min + 0.4*Lx, 5);
z_seeds = linspace(z_min + 0.1*Lz, z_min + 0.9*Lz, 3);

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

% Plot domain
plot([x(1) x(end) x(end) x(1) x(1)], ...
     [z(1) z(1)   z(end) z(end) z(1)], 'k-');

% Mark where turbine is
plot(x_turb,z_hub,'rx','Markersize',10,'Linewidth',2);

% Plot particle trajectories
for p = 1:numP
    plot(xp_all(p,:), zp_all(p,:), 'LineWidth', 1.2);
end

xlabel('X Position (m)');
ylabel('Z Position (m)');
title('Particle Trajectories');
xlim([x(1),x(end)]);
ylim([z(1),z(end)]);
grid on;
hold off;

