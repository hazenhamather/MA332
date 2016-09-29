function [xstar, fxstar, nitr, status] = Secant(f,a,b,epsilon,maxitr,loud)
%function [xstar, fxstar, nitr, status] = Secant(f,a,b,epsilon,maxitr,loud)
%This function uses the method of secants to search for the root of a
%function
%inputs
%f - function handle or filename
%a - left end-point of the interval
%b - right end-point of the interval
%epsilon - termination criteria for the dependent variable
%maxitr - a bound on the number of iterations
%loud - a flag indicating whether or not the algorithm should produce
%output on each iteratoin, with 0 indicating no output and nonzero
%indicating output
%Outputs
%xstar - the computed value of the root
%fxstar - the value of the function at the computed root
%nitr - the number of iterations used to compute the root
%status - the status of the search, with 0 indicating success and 1
%signaling failure to satisfy the convergence of tolerances

xstar = -Inf;
fxstar = -Inf;
status = 1;
nitr = 0;

if epsilon <= 0
    fprintf(2,'Please provide a positive epsilon value please.\n');
    return;
end

if abs(f(a)) < epsilon
    xstar = a;
    fxstar = f(a);
    status = 0;
    return;
end

if abs(f(b)) < epsilon
    xstar = b;
    fxstar = f(b);
    status = 0;
    return;
end

if a >= b
    fprintf(2,'Please make b larger than a.\n');
    return;
end

if abs(f(a)) < abs(f(b))
    xbest = abs(f(a));
else
    xbest = abs(f(b));
end

xstar = a-f(a)*(b-a)/(f(b)-f(a));

while abs(f(xstar)) > epsilon && nitr <= maxitr
   if f(b) == f(a)
       status = 1;
       fprintf(2,'Something went wrong.\n');
       return;
   end
   xstar = a-f(a)*(b-a)/(f(b)-f(a));
   fxstar = f(xstar);
   a = b;
   b = xstar;
   nitr = nitr + 1;
   if abs(f(a)) < abs(f(b)) && abs(f(a)) < xbest
        xbest = abs(f(a));
        xbestPrint = f(a);
    elseif abs(f(b)) < abs(f(a)) && abs(f(b)) < xbest
        xbest = abs(f(b));
        xbestPrint = f(b);
    end
   if loud
       fprintf('iter: %d; x_k = %f; best function value is %f\n',nitr,b,xbestPrint);
   end
    
end

if nitr > maxitr
    status = 1;
    fprintf(2,'Too many iterations.\n');
    return;
end

status = 0;