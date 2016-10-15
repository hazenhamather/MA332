function [L,status] = CholeskyInner(A)
%function [L,status] = CholeskyInner(A)
%This function computes a Cholesky Inner factorization where L*L'=A
%Input
%A-input matrix
%Outputs
%L-Lower triangular matrix
%status-states whether or not the process was successful; 0 on success, 1
%if not

m = size(A,1); %Getting number of rows
n = size(A,2); %Getting number of columns
status = 1; %Setting initial status to failure
L = zeros(m,n); %Initializing the Lower Triangular matrix
if m~=n %Making sure the input matrix is square
    fprintf(2,'Matrix must be square\n');
    return;
end
positivedefinite = all(eig(A) > 0); %Getting positive definite "boolean"
if m==1 && n==1 && A == 0 %Checking for zero
    L = 0;
    status = 0;
    return;
end
if ~positivedefinite %Making sure we have a positive definite input matrix
    fprintf(2,'Matrix must be positive definite\n');
    return;
end
if m == 1 && n == 1 %Checking for a 1x1 input matrix
    L = sqrt(A);
    status = 0;
    return;
end
for j = 1:m %Starting to move through the matrix, follow provided notes on 
            %Cholesky
    Ltranspose = ctranspose(L(j,1:j)); 
    temp = A(j,j)-L(j,1:j)*Ltranspose;
    if temp > 0
        L(j,j) = sqrt(temp);
        for i = j+1:m
            Ltrans = ctranspose(L(i,1:j));
            L(i,j) = A(j,i) - L(j,1:j)*Ltrans;
            L(i,j) = L(i,j)/L(j,j);
        end
    else
        return;
    end
end
status = 0;