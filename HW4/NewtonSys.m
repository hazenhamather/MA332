function [x, numitr, status] = NewtonSys(f,Jf,x0,epsilon,maxiter)
%function [x, numitr, status] = NewtonSys(f,Jf,x0,epsilon,maxiter)
%This function approximates the solution to a function by way of Newton's
%method through the Jacobian matrix. We iterate until either the conditions
%of epsilon or maxiter are met
%Inputs
%f-function to evaluate
%Jf - Jacobian matrix
%x0 - initial guess
%epsilon - user specified stopping criteria
%maxiter - the maximum number of iterations we will perform
%Outputs
%x - solution vector
%numitr - the number of iterations we perform
%status - says whether or not the function was a success, 1 on failure and
%0 on success

status = 1; %Setting initial status to failure
numitr = 0; %No iterations yet
x = x0; %Setting initial root to be infinity
[rows,~] = size(f);
% conditionMet = 0;
while 1 %Always in the loop until a stopping criterion is satisfied
    if rank(Jf(x)*f(x)) ~= rows %Checking for singularity
        fprintf(2,'Matrix is singular\n');
        return;
    end
    if all(abs(f(x)) < epsilon) %Checking for f(x) being less than epsilon
        status = 0;
        return;
    end
    if numitr >= maxiter %Checking if we have exceeded our max itr limit
        fprintf(2,'Maximum number of iterations reached\n');
        return;
    end
    deltax = -Jf(x)\f(x); %Get step size
    x = x + deltax; %Apply step size for the next iteration
    numitr = numitr + 1; %Update the iteration count
end