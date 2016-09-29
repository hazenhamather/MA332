function [v] = velocity(theta,m,alpha,w,d)
%function [v] = velocity(theta,m,alpha,w,d)
%This function computes the required velocity, v, to hit a target down
%range with a few input conditions
%inputs
%theta - initial launch angle
%m - mass of the projectile
%alpha - the drag coefficient
%w - wind velocity
%d - distance to target (in meters)
%outputs
%v - the initial velocity required in order to land the target downrange
v = fzero(@(x) projectile_landing(theta,x,m,alpha,w) - d, 3);
end