function [xstar, fxstar, nitr, status] = Newton(f,fp,x0,epsilon,maxitr,loud)
%function [xstar, fxstar, nitr, status] = Newton(f,fp,x0,epsilon,maxitr,loud)
%This function attempts to compute the root of a function, f, by way of
%Newton's method
%Inputs
%f - function handle or file for f(x)
%fp - a function handle or file for f'(x)
%x0 - the initial guess at the root
%epsilon - the termination criteria for the dependent variable
%maxitr - a bound on the number of iterations
%loud - a flag inidcating whether or not the alrogithm should produce
%output on each iteration, with 0 indicating no output and nonzero
%indicating output
%Outputs
%xstar - the computed value of the root
%fxstar - the function value at the computed root
%nitr - the number of iterations used to compute the root
%status - the status of the search, with 0 indicating success and 1
%signaling failure to satisfy the convergence tolerances

xstar = x0;
fxstar = f(xstar);
status = 1;
nitr = 0;

if epsilon <= 0
    fprintf(2,'Please supply a termination criteria greater than zero\n');
    return;
end

if abs(f(x0)) < epsilon
    status = 0;
    xstar = x0;
    fxstar = f(xstar);
    return;
end

while abs(f(xstar)) > epsilon && nitr <= maxitr
    if fp(xstar) ~= 0
        xstar = xstar - (f(xstar)/fp(xstar));
        if loud
            fprintf('Iter %d: Approximate solution: %f\n',nitr,f(xstar));
        end
%         xstar = xstar - (f(xstar)/fp(xstar));
    else
        fprintf(2,'The derivative was zero therefore I cannot continue\n');
        return;
    end
    nitr = nitr + 1;
end

if nitr > maxitr
    fprintf(2,'Too many iterations\n');
    status = 1;
    fxstar = f(xstar);
    return;
end

status = 0;
fxstar = f(xstar);
    