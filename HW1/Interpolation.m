function [xstar, fxstar, nitr, status] = Interpolation(f,a,b,epsilon,maxitr,loud);
%function [xstar, fxstar, nitr, status] = Interpolation(f,a,b,epsilon,maxitr,loud);
%This function applies Regula Falsi to a given function to identify a root
%inputs
%f - function filename or handle
%a - left end-point of the interval
%b - right end-point of the interval
%epsilon - terminatoin criteria for the dependent variable
%maxitr - a bound on the number of iterations
%loud - a flag indicating whether or not the algorithm should produce
%output on each iteration, with 0 indicating no output and nonzero
%indicating output
%outputs
%xstar - the computed value of the root
%fxstar - the value of the function evaluated at the computed root
%nitr - the number of iterations used to compute the root
%status - the status of the search, with 0 indicating success and 1
%indicating failure to satisfy the convergence tolerance
xstar = Inf;
fxstar = Inf;
status = 1;
nitr = 0;

if a > b
    fprintf(2, 'Sign error. Please make a larger than be and differ in sign\n');
    return;
end

if f(a) == 0 || abs(f(a)) < epsilon
    xstar = a;
    fxstar = f(xstar);
    status = 0;
    return;
end

if f(b) == 0 || abs(f(b)) < epsilon
    xstar = b;
    fxstar = f(xstar);
    status = 0;
    return;
end

if epsilon <= 0
    status = 1;
    return;
end

if abs(f(a)) < abs(f(b))
    oldbest = abs(f(a));
else
    oldbest = abs(f(b));
end

if loud
    fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n',nitr,a,b,oldbest);
end
xstar = b;
while abs(f(xstar)) > epsilon && nitr < maxitr
    xstar = a-f(a)*(a-b)/(f(a)-f(b));
    if sign(f(xstar)) == 0
        fxstar = f(xstar);
        status = 0;
        return;
    elseif sign(f(xstar)) == sign(f(a))
        a = xstar;
    elseif sign(f(xstar)) == sign(f(b))
        b = xstar;
    else
        status = 1;
        return;
    end
    
    nitr = nitr + 1;
 
    if loud
        if abs(f(a)) < abs(f(b)) && abs(f(a)) < oldbest
            newbest = abs(f(a));
            fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n',nitr,a,b,newbest);
            oldbest = newbest;
        elseif abs(f(b)) < abs(f(a)) && abs(f(b)) < oldbest
            newbest = abs(f(b));
            fprintf('Iteration %d: interval is [%f, %f]; best function value is %f\n',nitr,a,b,newbest);
            oldbest = newbest;
        end
        
    end

end

if nitr >= maxitr
    status = 1;
    nitr = nitr -1;
    fprintf(2,'Too many iterations\n');
    return;
end

nitr = nitr -1;
fxstar = f(xstar);
status = 0;