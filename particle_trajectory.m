% What this script does:
% given a velocity field u,w, it follows one particle through the flow in
% time and records its path

% It takes inputs:
% x0,z0- Initial particle position
% x,z- grids where u and w are defined
% u,w- Velocity fields
% dt- time step
% tmax- total simulation time
% U0- background horizontal velocity
% Lx,Lz- Boundary conditions

% It should return row vectors:
% xp- the x-position of the particle at each time step
% zp- the z position of the particle at each time step
% t- time vector

function[xp,zp,t]= particle_trajectory(x0,z0,x,z,u,w,U0,dt,tmax,Lx,Lz)
    % Number of time steps
    nt= round(tmax/dt)+1;
    
    % Set up initial values/memories
    xp = zeros(1, nt);
    zp = zeros(1, nt);
    t = linspace(0, tmax, nt);
    
    xp(1) = x0;
    zp(1) = z0;
    
    for n = 1:nt-1
        % Use interp_velocity to interpolate the velocity at the current position
        [up,wp] = interp_velocity(xp(n),zp(n),x,z,u,w);
    
        % Update particle position using differentiated versions of the eqns in Q6
        xp_new = xp(n) + (up + U0) * dt;
        zp_new = zp(n) + wp * dt;
    
        % Set values at upper/lower bounds
        if xp_new<0
            xp_new=xp_new+Lx;
        elseif xp_new>=Lx
            xp_new=xp_new-Lx;
        end
        
        if zp_new<0
            zp_new=0;
        elseif zp_new>=Lz
            zp_new=Lz;
        end
   
        xp(n+1) = xp_new;
        zp(n+1) = zp_new;
    end
end

