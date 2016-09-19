function [xstar, fxstar, nitr, status] = Bisection(f,a,b,epsilon,delta,maxitr,loud)
%function [xstar, fxstar, niter, status] = Bisection(f,a,b,epsilon,delta,maxitr,loud)
%This function attempts to find the root of a function using the bisection
%method. 
%Inputs
%f - the function to be evaluated
%a - the left end-point of the interval
%b - the right end-point of the interval
%epsilon - the termination tolerance for the independent variable
%delta - the termination tolerance for the dependent variable
%maxitr - a bound on the number of iterations
%loud - a flag indicating whether or not the algorithmm should produce
%output on each iteration, with 0 indicating no output and nonzero
%indicating output
%
%Outputs
%xstar - the computed value of the root
%fxstar - the value of the function at the computed root
%niter - the number of iterations used to compute the root
%status - the status of the search, with 0 indicating success and 1
%signaling failure to satisfy one of the convergence tolerances

xstar = -Inf;
fxstar = -Inf;
status = 1;
nitr = 0;

if f(a)*f(b) > 0
    fprintf(2, 'I need a bracket\n');
    return;
end

if delta <= 0
    fprintf(2, 'I need a positive stopping criteria.\n');
    return;
end

if epsilon <= 0
    fprintf(2, 'I need a positive stopping criteria.\n');
    return;
end

if abs(f(a)) < epsilon
    status = 0;
    xstar = a;
    fxstar = f(xstar);
    return;
end

if abs(f(b)) < epsilon
    status = 0;
    xstar = b;
    fxstar = f(xstar);
    return;
end

if abs(f(a)) < abs(f(b))
    xbest = abs(f(a));
else
    xbest = abs(f(b));
end

% if loud
%     fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n', nitr,a,b,xbest);
% end

x = (a+b)/2;

while (abs(f(x)) > epsilon || abs(b-a) > delta) && nitr <= maxitr
    if loud
        fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n', nitr,a,b,xbest);
    end
    x = (a+b)/2;
    xstar = x;
    fxstar = f(xstar);
    if fxstar == 0
        status = 0;
        return;
    end
    
    if f(a) * fxstar < 0
%         a = a;
        b = xstar;
    else
        a = xstar;
%         b = b;
    end
    
    if abs(f(a)) < abs(f(b)) && abs(f(a)) < xbest
        xbest = abs(f(a));
    elseif abs(f(b)) < abs(f(a)) && abs(f(b)) < xbest
        xbest = abs(f(b));
    end
  
    nitr = nitr + 1;
    
%     if loud
%         fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n', nitr,a,b,xbest);
%     end
    
end

if nitr > maxitr
    fprintf(2,'Too many iterations.\n');
    return;
end

status = 0;

% 
% if ~exist('loud')
%     loud = 0;
% end
% 
% if f(a) == 0
%     xstar = a;
%     fxstar = f(xstar);
%     status = 0;
%     return;
% end
% 
% if f(b) == 0
%     xstar = b;
%     fxstar = f(xstar);
%     status = 0;
%     return;
% end
% if loud ~= 0
%     fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n', nitr,a,b,b);
% end
% 
% c = (a+b)/2;
% xstar = c;
% fxstar = f(xstar);
% % if nitr > maxitr
% %     xstar = b;
% %     fxstar = f(xstar);
% %     return;
% % end
% while abs(b-a) > epsilon && abs(fxstar) > delta
%     c = (a+b)/2;
%     xstar = c;
%     fxstar = f(xstar);
%     
%     if fxstar == 0
%         status = 0;
%         xstar = c;
%         nitr = nitr + 1;
%         return;
%     end
%     if f(a)*f(xstar) < 0
%         a = a;
%         b = c;
%     else
%         a = c;
%         b = b;
%     end
%     nitr = nitr + 1;
%     if loud ~= 0
%         fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n', nitr,a,b,fxstar);
%     end
% end
% 
% status = 0;

    
    