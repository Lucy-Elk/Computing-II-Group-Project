%% 
% Wind Turbine Flow Solver - Tasks 4 & 5
% Spectral method in x, finite differences in z

clear; clc;

% Setup parameters
Lx = 1000;          
Lz = 800;
Nx = 256;          
Nz = 128;           
U0 = 20;          

dx = Lx / Nx;
dz = Lz / (Nz + 1);
x = linspace(0, Lx - dx, Nx);
z = linspace(0, Lz, Nz + 2);  % including ghost points


% Define forcing - using gaussian to model turbine
x_0 = Lx/4;     
z_0 = 150;       
sigma_x = 50;           
sigma_z = 40; 

% Choose a large A to amplify the effect
A = 5;  

[X, Z] = meshgrid(x, z);
fx = A * exp(-((X - x_0).^2 / (2 * sigma_x^2) + ...
                       (Z - z_0).^2 / (2 * sigma_z^2)));

% FFT in x-direction
fx_hat = fft(fx, [], 2) / Nx;  
w_hat = zeros(Nz + 2, Nx);

% pis because its periodic
k_vec = 2 * pi * [0:Nx/2, -Nx/2+1:-1] / Lx;

% Main loop - solve for each wavenumber
for k = 1:Nx
    bk = k_vec(k);
    
    % Build matrix A - same as derived in task 3
    A = zeros(Nz + 2, Nz + 2);
    
    A(1, 1) = 1;  % bottom BC: w=0
    
    % interior poits
    for j = 2:Nz+1
        A(j, j-1) = 1;
        A(j, j) = -(2 + dz^2 * bk^2);
        A(j, j+1) = 1;
    end
    
    % top BC: dw/dz = 0
    A(Nz+2, Nz) = -1;
    A(Nz+2, Nz+1) = 0;
    A(Nz+2, Nz+2) = 1;
    
    % RHS vector
    rhs = zeros(Nz + 2, 1);
    %b(1) = 0;
    
    for j = 2:Nz+1
        rhs(j) = dz * (fx_hat(j+1, k) - fx_hat(j-1, k)) / (2 * U0);
    end
    
    %b(Nz+2) = 0;
    
    % solve system
    w_hat(:, k) = A \ rhs;
    
    if k==1
        w_hat(:,k)=0;
    end
end

% inverse FFT to get w(x,z)
w = real(ifft(w_hat * Nx, [], 2));

% Recover u from continuity: du/dx = -dw/dz
dwdz = zeros(Nz + 2, Nx);
for j = 2:Nz+1
    dwdz(j, :) = (w(j+1, :) - w(j-1, :)) / (2 * dz);
end
dwdz(1, :) = (w(2, :) - w(1, :)) / dz;
dwdz(Nz+2, :) = (w(Nz+2, :) - w(Nz+1, :)) / dz;

% FFT then integrate by dividing by ik
dwdz_hat = fft(dwdz, [], 2) / Nx;
u_hat = zeros(Nz + 2, Nx);
for k = 1:Nx
    if abs(k_vec(k)) > 1e-10  
        u_hat(:, k) = 1i * dwdz_hat(:, k) / k_vec(k);
    else
        u_hat(:, k) = 0;  
    end
end

u = real(ifft(u_hat * Nx, [], 2));
u=u-u(:,1);

% Plots
figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
pcolor(X, Z, w);
shading interp;
colorbar;
xlabel('x (m)');
ylabel('z (m)');
title('Vertical velocity w(x,z)');
colormap(jet);

subplot(2, 3, 2);
pcolor(X, Z, u);
shading interp;
colorbar;
xlabel('x (m)');
ylabel('z (m)');
title('Horizontal velocity u(x,z)');
colormap(jet);

subplot(2, 3, 3);
contourf(X, Z, w, 20);
colorbar;
xlabel('x (m)');
ylabel('z (m)');
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
xlabel('x (m)');
ylabel('w (m/s)');
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
xlabel('x (m)');
ylabel('u (m/s)');
title('u at different heights');
legend('Location', 'best');
grid on;

subplot(2, 3, 6);
pcolor(X, Z, fx);
shading interp;
colorbar;
xlabel('x (m)');
ylabel('z (m)');
title('Forcing function f_x(x,z)');
colormap(jet);

sgtitle('Wind Turbine Flow Simulation Results');

%% save
save('wind_turbine_results.mat', 'x', 'z', 'w', 'u', 'fx', 'X', 'Z','U0','Lx','Lz','x_0','z_0');
fprintf('Done!\n');

