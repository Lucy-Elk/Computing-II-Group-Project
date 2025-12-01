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
    xp = min(max(xp, x(1)), x(end));
    zp = min(max(zp, z(1)), z(end));

    % making everything the form needed to use interp2
    [Xgrid,Zgrid]= meshgrid(x,z);
    
    % Use interp2 to estimate the velocity between nodes to mkae particles
    % move continuously
    up= interp2(Xgrid,Zgrid,u,xp,zp,'linear',0);
    wp= interp2(Xgrid,Zgrid,w,xp,zp,'linear',0);
end
