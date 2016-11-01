function [x,fx,y,s,status,niter]=interior_point_quadratic(Q,c,A,b,x0,y0,s0,tol,maxiter,loud)
%[x,fx,y,s,status,niter]=interior_point_linear(Q,c,A,b,x0,y0,s0,tol,maxiter,loud)
%Solve min_x Qx+c'x 
%subject to
%      Ax=b
%      x>=0
%using an interior point method. 
%
%INPUT
%
%Q,c,A,b -- from the optimization problem
%x0    -- initial x (default to all 1's)
%y0    -- inital guess at equality Lagrange multipliers (default to all
%1's)
%s0    -- initial guess at inequality Lagrange multipliers (default to all
%1's)
%tol   -- Iterate until solution found with s'*x <= tol (default 1e-6)
%maxiter -- maximum number of iterations to perform (default 1e3)
%loud  -- whether or not to produce output.  loud==0 gives no output,
%loud==1 gives modest output at each iteration, loud == 2 gives lots of
%output at each iteration. 
%
%OUTPUT
%
%x     -- approximate minimizer
%fx    -- c'*x
%y     -- equality Lagrange multiplier
%s     -- inequality Lagrange multiplier
%status -- 0 means success, nonzero means failure
%niter -- number of iterations required for convergence. 

n = size(c,1);
e = ones(n,1);

if ~exist('x0') || isempty(x0)
    x0 = e;
end
if ~exist('y0') || isempty(y0)
    y0 = ones(size(A,1),1);
end
if ~exist('s0') || isempty(s0)
    s0 = e;
end
if~exist('tol') || isempty(tol)
    tol = 1e-6;
end
if ~exist('maxiter') || isempty(maxiter)
    maxiter = 1e6;
end
if ~exist('loud') || isempty(loud)
    loud = 0;
end


%initialize
x=x0;
fx=0.5*x'*Q*x+c'*x;
y=y0;
s=s0;
status=1;
niter=Inf;

dig=3; %number of digits to use when producing output. 


%Make sure everything is legit:
%c must be a column vector
assert(size(c,2)==1,'c must be a column vector');
n=size(c,1);  %n is the number of design variables. 

%A must be something by n
assert(size(A,2)==n,'A must have number of columns equal to number of design variables, in this case %d',n);
m=size(A,1); %m is the number of equality constratints. 

%b must be a column vector of shape compatible with Ax=b. 
assert(size(b,1)==m && size(b,2)==1,'b must be a column vector of shape compatible with Ax=b');

%x0 must be non-negative
assert(all(x0>=0), 'x0 must be non-negative');

%x0 must have shape compatbile with A
assert(size(x0,1)==n && size(x0,2)==1,'x0 must be a column vector with shape compatible with A');

%y0 must be m x 1
assert(size(y0,1)==m && size(y0,2)==1,'y0 must have as many rows as equality contraints and exactly 1 column.');

%s0 must be n x 1
assert(size(s0,1)==n && size(s0,2)==1,'s0 must have as many rows as design variables and exactly 1 column.');

%s0 must be non-negative
assert(all(s0>=0),'s0 must be non-negative');

%tol must be postibe
assert(tol>0,'tol must be positive'); 

%maxiter must be positive
assert(maxiter>0,'maxiter must be positive.')





%Make a diagonal matrix out of the vector. 
S=@(s) diag(s);
X=@(x) diag(x);

%The Lagrangian is
%L(x,y,s) = c' x - y' (Ax-b) -s' x

%The KKT conditions are then
%  c-A'y-s=0
%  Ax-b   =0
%  s'x    =0 , s>=0, x>=0

%These are relaxed to 

%Ax-b=0
%A'y+s-c=0
%Xs=mu e
%at each step.

%Define a function 
%F(x,y,s)=[Ax-b; A'y+s-c;Xs-mu e]
%Then 
%JF=[A, 0, 0; 0 A' I; S 0 X]
%We can use Newton's method to solve at each step. 




F=@(x,y,s,mu) [A*x-b; -Q*x+A'*y+s-c; X(x)*s-mu*e];

JF=@(x,y,s) [A,              zeros(m,m)          zeros(m,n);...
             -Q,     A',                 eye(n,n); 
             S(s),           zeros(n,m),         X(x)];

mu=s'*x;

N=1;
while s'*x>tol && N<maxiter
    %Take one step of newton's method
    %the solution to JF*Delta=-F will contain the updates [deltax, deltay,
    %deltas.
    Delta=-JF(x,y,s)\F(x,y,s,mu);
    deltax=Delta(1:n);
    deltay=Delta(n+1:n+m);
    deltas=Delta(n+m+1:end);
    
    %If JF was singular, it is possible that we are all screwed up now.
    %Check to make sure this didn't happen.
    
    if any(isnan(Delta)) || any(isinf(Delta))
        return;
    end
    
    
    %going to set x=x+alphax*deltax.  Make sure we don't force x or s to be
    %infeasible. 
    
    %If all of deltax, deltas are non-negative, then go ahead and use full
    %newton step
    if all(deltax>=0) && all(deltas>=0)
        alpha=1;
    else
        alphax=min(x(deltax<0)./abs(deltax(deltax<0)))*0.95;
        alphas=min(s(deltas<0)./abs(deltas(deltas<0)))*0.95;
        alpha=min([alphax,alphas]);
    end
    
    
    %A crazy time amount of output for debugging what is up.
    if loud >= 2
        fprintf('\n\n\n');
        fprintf('k: %d x_k: %s y_k: %s s_k: %s\n',N,mat2str(x,dig), mat2str(y,dig), mat2str(s,dig));
        fprintf('      mu = %s\n',mat2str(mu,dig));
        fprintf('      JF = %s\n',mat2str(JF(x,y,s),dig));
        fprintf('      -F = %s\n',mat2str(-F(x,y,s,mu),dig));
        fprintf('      [deltax,deltay,deltas]=%s\n',mat2str(Delta,dig));
        fprintf('      alphax=%s , alphas=%s\n',mat2str(alphax,dig),mat2str(alphas,dig));
        fprintf('      x_k+1=%s, y_k+1=%s, s_k+1=%s\n',mat2str(x+alpha*deltax,dig), mat2str(y+alpha*deltay,dig), mat2str(s+alpha*deltas,dig));
    end

        
    %A moderate amount of output for monitoring progress
    if loud >= 1
        %This is a nice way to print out formatted vectors; use mat2str to
        %make a nice string then use %s to just print a string.
        fprintf('Iteration: %d x=%s y=%s s=%s\tx''s=%s ||Ax-b||=%s ||c-A''y-s||=%s\n',...
                             N,mat2str(x,dig),mat2str(y,dig),mat2str(s,dig),...
                             mat2str(x'*s,dig),mat2str(max(abs(A*x-b)),dig),mat2str(max(abs(c-A'*y-s)),dig));
    end
    %Now update 
    x=x+alpha*deltax;
    y=y+alpha*deltay;
    s=s+alpha*deltas;
    
    fx=0.5*x'*Q*x+c'*x;
    
    %Turn down mu
    mu=s'*x/N^2;
    
    N=N+1;
end
    
if N<maxiter
    status=0;
    niter=N;
end


