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

%Setting some variables
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

%Checking for a correct initial guess
if norm(gf(xstar)) < gradTol
    return;
end
%Making sure ireBound was not zero or negative
if numitr >= itrBound
    fprintf(2,'Exceeded number of allowed iterations\n');
    return;
end

%Entering Steepest Descent
if dparam == 1
    shouldContinue = true; %Stopping variable
    xkOld = x0;
    while shouldContinue
        dk = -gf(xkOld);
        if dk == 0
            return;
        end
        switch lparam %checking for which lparam we were given
            case 1 %constant alpha
                alpha = 1;
            case 2 %interpolating polynomial for alpha
                alpha = (-3*f(xkOld)+4*f(xkOld+dk)-f(xkOld+2*dk))/(-2*f(xkOld)+...
                4*f(xkOld+dk)-2*f(xkOld+2*dk));
            case 3 %bisection method for alpha if needed
                if gf(xkOld)*gf(xkOld + dk) > 0
                    alpha = 1;
                else
                    [alpha,~,~] = alphaBisection(gf,dk,xkOld,0,1);
                end
            otherwise %This will trip if the user provides an incorrect lparam
                fprintf(2,'lparam must be between 1 and 3 inclusive\n');
                return;
        end
        xkNew = xkOld + alpha * dk; %getting next x_k
        xtempOld = xkOld; %setting temporary variable for old x_k
        xkOld = xkNew;
        numitr = numitr + 1;
        switch loud %does the user want output?
            case {1,2}
                if lparam == 1
                    if numitr - 1 == 0
                        fprintf('Optimization method:\n');
                        fprintf('Search Direction: Steepest Descent\n');
                        fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                    end
                else
                    fprintf('Optimization method:\n');
                    fprintf('Search Direction: Steepest Descent\n');
                    fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                end
        end
        switch loud  %more output            
            case 2
                fprintf('Iteration %d\n',numitr-1);
                fprintf('x_k = (%f)\n',xtempOld);
                fprintf('d = (%f)\n',dk);
                fprintf('alpha = %f\n',alpha);
                fprintf('x_k+1 = (%f)\n',xkNew);
                fprintf('norm(grad f)=%f, xdiff=%f\n\n\n',norm(gf(xkNew)),abs(xkNew-xtempOld));
        end
        if norm(gf(xkNew)) < gradTol || abs(xkNew - xtempOld) < xTol
            %Checking tolerances to see if we satisfied any criteria ^
            status = 0;
            xstar = xkNew;
            fxstar = f(xstar);
            gradientCriteriaMet = true; %Flag for final output
            break;
        end
        if numitr > itrBound %did we reach our iteration limit?
            xstar = xkNew;
            fxstar = f(xstar);
            iterationLimitReached = true; %Flag for final output
            break;
        end
    end
    if loud == 1 || loud == 2 %Final output if specified by user
        fprintf('Optimization Terminated\n');
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
    
    %The above comments in the Steepest Descent alrogithm apply to the next
    %two. In retrospect, after writing these comments, I realized that I
    %could have written this code to be about 1/3 as long as it currently
    %is. I should have put all of the initial print statements at the top
    %of the file before entering the methods. This would have allowed me to
    %execute methods, store every result into a vector, and print output
    %from that vector. If memory or readability was an immediate issue, I
    %would change the code as I proposed but as of right now, I have plenty of
    %memory and seeing as this code is approximately 13kB in size, I doubt
    %many people are hurting for it.
    
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
            case {1,2}
                if lparam == 1
                    if numitr - 1 == 0
                        fprintf('Optimization method:\n');
                        fprintf('Search Direction: Steepest Descent\n');
                        fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                    end
                else
                    fprintf('Optimization method:\n');
                    fprintf('Search Direction: Steepest Descent\n');
                    fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                end
        end
        switch loud              
            case 2
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
                    break;
                end
    end
    if loud == 1 || loud == 2
        fprintf('Optimization Terminated\n');
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
            case {1,2}
                if lparam == 1
                    if numitr - 1 == 0
                        fprintf('Optimization method:\n');
                        fprintf('Search Direction: Steepest Descent\n');
                        fprintf('Step Length Selection: Constant alpha = 1\n\n\n');
                    end
                else
                    fprintf('Optimization method:\n');
                    fprintf('Search Direction: Steepest Descent\n');
                    fprintf('Step Length Selection: alpha = %f\n\n\n',alpha);
                end
        end
        switch loud              
            case 2
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
            break;
        end
    end
    if loud == 1 || loud == 2
        fprintf('Optimization Terminated\n');
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