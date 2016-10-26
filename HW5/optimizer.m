function [xstar,fxstar,status] = optimizer(f,gf,Hf,x0,gradTol,xTol,itrBound,dparam,lparam,loud)
%function [xstar,fxstar,status] = optimizer(f,gf,Hf,x0,gradTol,xTol,itrBound,dparam,lparam,loud)
%This function optimizes a given function, f, with parameters and has the
%potential to output details of the process along the way. This solution
%utilizes a companion bisection solver function
%Inputs
%f - the function to be minimized
%gf - the gradient of f
%Hf - the hession of the f
%x0 - an initial guess or iterate
%gradTol - a termination tolerance, iteration should cease successfully if
    %gf(xk) less than this value
% xTol - a termination tolerance, iteration should stop successfully if the
%   abs(xk-(xk-1)) is less than this value
%itrBound - an iteration limit, optimization will be unsuccussful if this
%   value is met
%dparam - a direction choice parameter
%lparam - a line search parameter
%loud - a flag indicating whether or not to print diagnostic information
%Outputs
%xstar - the last computed xk, best estimate of the optimization
%fxstar - the function value at xstar
%status - a status variable that will be 0 if the algorithm terminates with
%   either the gradient or iteration tolerance satisfied, will be nonzero
%   otherwise
status = 1;
xstar = x0;
fxstar = f(xstar);
numitr = 0;
gradientCriteriaMet = false;
iterationLimitReached = false;
if ~exist('Hf') || isempty(Hf)
    Hf = 1;
end
if ~exist('gradTol') || isempty(gradTol) || gradTol <= 0
    gradTol = 1e-6;
end
if ~exist('xTol') || isempty(xTol) || xTol <= 0
    xTol = 1e-6;
end
if ~exist('itrBound') || isempty(itrBound) || itrBound <= 0
    itrBound = 1e2;
end
if ~exist('dparam') || isempty(dparam) || dparam <= 0
    dparam = 1;
end
if ~exist('lparam') || isempty(lparam) || lparam <= 0
    lparam = 2;
end
if ~exist('loud') || isempty(loud) || loud < 0
    loud = 0;
end
if norm(gf(xstar)) < gradTol
    return;
end
if numitr >= itrBound
    fprintf(2,'Exceeded number of allowed iterations\n');
    return;
end
if dparam == 1
    shouldContinue = true;
    xkOld = x0;
    while shouldContinue
        dk = -gf(xkOld);
        if dk == 0
            return;
        end
        switch lparam
            case 1
                alpha = 1;
            case 2
                alpha = (-3*f(xkOld)+4*f(xkOld+dk)-f(xkOld+2*dk))/(-2*f(xkOld)+...
                4*f(xkOld+dk)-2*f(xkOld+2*dk));
            case 3
                if gf(xkOld)*gf(xkOld + dk) > 0
                    alpha = 1;
                else
                    [alpha,~,~] = alphaBisection(gf,dk,xkOld,0,1);
                end
            otherwise
                fprintf(2,'lparam must be between 1 and 3 inclusive\n');
                return;
        end
        xkNew = xkOld + alpha * dk;
        xtempOld = xkOld;
        xkOld = xkNew;
        numitr = numitr + 1;
        switch loud
            case 1
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end               
            case 2
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end 
                fprintf('Iteration %d\n',numitr-1);
                fprintf('x_k = (%f)\n',xtempOld);
                fprintf('d = (%f)\n',dk);
                fprintf('alpha = %f\n',alpha);
                fprintf('x_k+1 = (%f)\n',xkNew);
                fprintf('norm(grad f)=%f, xdiff=%f\n\n\n',norm(gf(xkNew)),abs(xkNew-xtempOld));
        end
        if norm(gf(xkNew)) < gradTol || abs(xkNew - xtempOld) < xTol
            status = 0;
            xstar = xkNew;
            fxstar = f(xstar);
            gradientCriteriaMet = true;
            break;
        end
        if numitr > itrBound
            xstar = xkNew;
            fxstar = f(xstar);
            iterationLimitReached = true;
            return;
        end
    end
    if loud == 1 || loud == 2
        fprintf('Operation Terminated\n');
        fprintf('Total Iterations: %d\n', numitr);
        if status == 0
            fprintf('Status: Successful\n');
            if gradientCriteriaMet == true
                fprintf('Reason: Gradient norm criterion reached\n\n\n');
            end
        else
            fprintf('Status: Failed\n');
            if iterationLimitReached == true
                fprintf('Maximum number of iterations reached\n\n\n');
            end
        end
    end
%Newton Direction
elseif dparam == 2
    shouldContinue = true;
    xkOld = x0;
    while shouldContinue
        dk = -Hf(xkOld)\gf(xkOld);
        if dk == 0
            return;
        end
        switch lparam
            case 1
                alpha = 1;
            case 2
                alpha = (-3*f(xkOld)+4*f(xkOld+dk)-f(xkOld+2*dk))/(-2*f(xkOld)+...
                4*f(xkOld+dk)-2*f(xkOld+2*dk));
            case 3
                if gf(xkOld)*gf(xkOld + dk) > 0
                    alpha = 1;
                else
                    [alpha,~,~] = alphaBisection(gf,dk,xkOld,0,1);
                end
            otherwise
                fprintf(2,'lparam must be between 1 and 3 inclusive\n');
                return;
        end
        xkNew = xkOld + alpha * dk;
        xtempOld = xkOld;
        xkOld = xkNew;
        numitr = numitr + 1;
        switch loud
            case 1
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end               
            case 2
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end
        end
                
                fprintf('Iteration %d\n',numitr-1);
                fprintf('x_k = (%f)\n',xtempOld);
                fprintf('d = (%f)\n',dk);
                fprintf('alpha = %f\n',alpha);
                fprintf('x_k+1 = (%f)\n',xkNew);
                fprintf('norm(grad f)=%f, xdiff=%f\n\n\n',norm(gf(xkNew)),abs(xkNew-xtempOld));
                if norm(gf(xkNew)) < gradTol || abs(xkNew - xtempOld) < xTol
                    status = 0;
                    xstar = xkNew;
                    fxstar = f(xstar);
                    break;
                end
                if numitr > itrBound
                    xstar = xkNew;
                    fxstar = f(xstar);
                    break;
                end
    end
    if loud == 1 || loud == 2
        fprintf('Operation Terminated\n');
        fprintf('Total Iterations: %d\n', numitr);
        if status == 0
            fprintf('Status: Successful\n');
            if gradientCriteriaMet == true
                fprintf('Reason: Gradient norm criterion reached\n\n\n');
            end
        else
            fprintf('Status: Failed\n');
            if iterationLimitReached == true
                fprintf('Maximum number of iterations reached\n\n\n');
            end
        end
    end
%BFGS
elseif dparam == 3
    shouldContinue = true;
    xkOld = x0;
    while shouldContinue
        dk = -Hf\gf(xkOld);
        if dk == 0
            return;
        end
        switch lparam
            case 1
                alpha = 1;
            case 2
                alpha = (-3*f(xkOld)+4*f(xkOld+dk)-f(xkOld+2*dk))/(-2*f(xkOld)+...
                4*f(xkOld+dk)-2*f(xkOld+2*dk));
            case 3
                if gf(xkOld)*gf(xkOld + dk) > 0
                    alpha = 1;
                else
                    [alpha,~,~] = alphaBisection(gf,dk,xkOld,0,1);
                end
            otherwise
                fprintf(2,'lparam must be between 1 and 3 inclusive\n');
                return;
        end
        xkNew = xkOld + alpha * dk;
        xtempOld = xkOld;
        xkOld = xkNew;
        numitr = numitr + 1;
        deltagf = gf(xkNew) - gf(xtempOld);
        deltax = xkNew - xtempOld;
        Hf = Hf + ((deltagf*deltagf')/(deltagf'*deltax))-((Hf*(deltax*deltax')*Hf)/(deltax'*Hf*deltax));
        switch loud
            case 1
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end               
            case 2
                switch lparam
                    case 1
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                        end
                    case 2
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                    case 3
                        if numitr - 1 == 0
                            fprintf('Optimization method:\n');
                            fprintf('Search Direction: Steepest Descent\n');
                            fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                        end
                end 
                fprintf('Iteration %d\n',numitr-1);
                fprintf('x_k = (%f)\n',xtempOld);
                fprintf('d = (%f)\n',dk);
                fprintf('alpha = %f\n',alpha);
                fprintf('x_k+1 = (%f)\n',xkNew);
                fprintf('norm(grad f)=%f, xdiff=%f\n\n\n',norm(gf(xkNew)),abs(xkNew-xtempOld));
        end
        if norm(gf(xkNew)) < gradTol || abs(xkNew - xtempOld) < xTol
            status = 0;
            xstar = xkNew;
            fxstar = f(xstar);
            break;
        end
        if numitr > itrBound
            xstar = xkNew;
            fxstar = f(xstar);
            return;
        end
    end
    if loud == 1 || loud == 2
        fprintf('Operation Terminated\n');
        fprintf('Total Iterations: %d\n', numitr);
        if status == 0
            fprintf('Status: Successful\n');
            if gradientCriteriaMet == true
                fprintf('Reason: Gradient norm criterion reached\n\n\n');
            end
        else
            fprintf('Status: Failed\n');
            if iterationLimitReached == true
                fprintf('Maximum number of iterations reached\n\n\n');
            end
        end
    end
end
end