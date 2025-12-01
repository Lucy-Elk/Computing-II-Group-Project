% Wind Turbine Flow Solver - Tasks 4 & 5
% Spectral method in x, finite differences in z

clear; close all; clc;

%% Setup parameters
Lx = 10;          
H = 5;             
Nx = 128;          
Nz = 64;           
U0 = 1.0;          

dx = Lx / Nx;
dz = H / (Nz + 1);

x = linspace(0, Lx - dx, Nx);
z = linspace(0, H, Nz + 2);  % including ghost points

%% Define forcing - using gaussian to model turbine
x_turbine = Lx / 2;     
z_turbine = H / 3;       
sigma_x = 0.5;           
sigma_z = 0.8;           
amplitude = -2.0;        

[X, Z] = meshgrid(x, z);
fx = amplitude * exp(-((X - x_turbine).^2 / (2 * sigma_x^2) + ...
                       (Z - z_turbine).^2 / (2 * sigma_z^2)));

%% FFT in x-direction
fx_hat = fft(fx, [], 2) / Nx;  
w_hat = zeros(Nz + 2, Nx);

k = 2 * pi * [0:Nx/2, -Nx/2+1:-1] / Lx;

%% Main loop - solve for each wavenumber
for m = 1:Nx
    km = k(m);
    
    % Build matrix A - same as derived in task 3
    A = zeros(Nz + 2, Nz + 2);
    
    A(1, 1) = 1;  % bottom BC: w=0
    
    % interior poits
    for j = 2:Nz+1
        A(j, j-1) = 1;
        A(j, j) = -(2 + dz^2 * km^2);
        A(j, j+1) = 1;
    end
    
    % top BC: dw/dz = 0
    A(Nz+2, Nz) = -1;
    A(Nz+2, Nz+1) = 0;
    A(Nz+2, Nz+2) = 1;
    
    % RHS vector
    b = zeros(Nz + 2, 1);
    %b(1) = 0;
    
    for j = 2:Nz+1
        b(j) = dz * (fx_hat(j-1, m) - fx_hat(j+1, m)) / (2 * U0);
    end
    
    %b(Nz+2) = 0;
    
    % solve system
    w_hat(:, m) = A \ b;
end

%% inverse FFT to get w(x,z)
w = real(ifft(w_hat * Nx, [], 2));

%% Recover u from continuity: du/dx = -dw/dz
dwdz = zeros(Nz + 2, Nx);
for j = 2:Nz+1
    dwdz(j, :) = (w(j+1, :) - w(j-1, :)) / (2 * dz);
end
dwdz(1, :) = (w(2, :) - w(1, :)) / dz;
dwdz(Nz+2, :) = (w(Nz+2, :) - w(Nz+1, :)) / dz;

% FFT then integrate by dividing by ik
dwdz_hat = fft(dwdz, [], 2) / Nx;
u_hat = zeros(Nz + 2, Nx);
for m = 1:Nx
    if abs(k(m)) > 1e-10  
        u_hat(:, m) = 1i * dwdz_hat(:, m) / k(m);
    else
        u_hat(:, m) = 0;  
    end
end

u = real(ifft(u_hat * Nx, [], 2));

%% Plots
figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
pcolor(X, Z, w);
shading interp;
colorbar;
xlabel('x');
ylabel('z');
title('Vertical velocity w(x,z)');
colormap(jet);

subplot(2, 3, 2);
pcolor(X, Z, u);
shading interp;
colorbar;
xlabel('x');
ylabel('z');
title('Horizontal velocity u(x,z)');
colormap(jet);

subplot(2, 3, 3);
contourf(X, Z, w, 20);
colorbar;
xlabel('x');
ylabel('z');
title('Contours of w(x,z)');

% plot w at a few different heights
subplot(2, 3, 4);
heights_idx = round([Nz/4, Nz/2, 3*Nz/4]) + 1;
hold on;
for idx = heights_idx
    plot(x, w(idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('z = %.2f', z(idx)));
end
hold off;
xlabel('x');
ylabel('w');
title('w at different heights');
legend('Location', 'best');
grid on;

subplot(2, 3, 5);
hold on;
for idx = heights_idx
    plot(x, u(idx, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('z = %.2f', z(idx)));
end
hold off;
xlabel('x');
ylabel('u');
title('u at different heights');
legend('Location', 'best');
grid on;

subplot(2, 3, 6);
pcolor(X, Z, fx);
shading interp;
colorbar;
xlabel('x');
ylabel('z');
title('Forcing function f_x(x,z)');
colormap(jet);

sgtitle('Wind Turbine Flow Simulation Results');

%% save
save('wind_turbine_results.mat', 'x', 'z', 'w', 'u', 'fx', 'X', 'Z','U0','Lx','H','x_turbine','z_turbine');
fprintf('Done!\n');

