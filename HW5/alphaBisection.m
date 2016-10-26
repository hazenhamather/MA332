function [alpha,nitr,status] = alphaBisection(gf,dk,x0,a,b,epsilon,delta,maxitr)
%function [xstar,fxstar,nitr,status] = alphaBisection(f,gf,a,b,epsilon,delta,maxitr,loud)
%This function will solve for alpha via bisection when given the function
%to be minimized along with its gradient.
%Inputs
%f - function to be evaluated
%gf - gradient of f
%a - left endpoint of the possible alpha values
%b - right endpoint of the possible alpha values
%epsilon - termination toerlance for the independent variable
%delta - termination tolerance for the dependent variable
%maxitr - a bound on the number of iterations
%loud - a flag indicating whether or not the algorithm should produce
%output on each iteration, with 0 indicating no output and nonzero
%indicating output
%Outputs
%xstar - the computed value of the root
%fxstar - the value of the function at the computed root
%nitr - the number of iterations used to compute the root
%status - the status of the search, with 0 indicating success and 1
%signaling failure to satisfy one of the convergence tolerances

status = 1;
nitr = 0;

if ~exist('epsilon') || isempty(epsilon)
    epsilon = 1e-12;
end
if ~exist('delta') || isempty(delta)
    delta = 1e-12;
end
if ~exist('maxitr') || isempty(maxitr)
    maxitr = 1000;
end

xkOld = x0;
fxstar = xkOld;
alpha = (a+b)/2;
while (abs(fxstar) > epsilon || abs(b-a) > delta) && nitr <= maxitr
    alpha = (a+b)/2;
    fxstar = gf(xkOld + alpha*dk)*dk;
    if fxstar == 0
        status = 0;
        return;
    end
    if gf(a) * fxstar < 0
        b = alpha;
    else
        a = alpha;
    end
    nitr = nitr + 1;
end
if nitr > maxitr
    return;
end
status = 0;
end