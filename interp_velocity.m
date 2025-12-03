% What this script does:
% Takes a particle's current position (xp,zp),
% and finds the value of the u(x,z) and w(x,z) at that exact point.
% Since the particles position has only previously been defined on a grid, 
% the function uses interpolation to estimate the velocity from the
% surrounding grid values.

% It takes inputs:
% (xp,zp)- the particle's current position
% (x,z) vectors of grid values 
% u- horizontal velocity 
% w- vertical velocity 

% It should return:
% up- interpolated horizontal velocity at (xp,zp)
% wp- interpolated vertical velocity at (xp,zp)

function[up,wp]=interp_velocity(xp,zp,x,z,u,w)
    dx= x(2)-x(1);
    Lx = x(end) - x(1) +dx;
    xp=x(1)+mod(xp-x(1),Lx);

    zp=min(max(zp,z(1)),z(end));
    % Use interp2 to estimate the velocity between nodes to mkae particles
    % move continuously
    up= interp2(x,z,u,xp,zp,'linear',0);
    wp= interp2(x,z,w,xp,zp,'linear',0);
end
