function [L,U,status] = lu_nopivot(A)
%function [L,U,p] = lu_nopivot(A)
%The purpose of this function is to compute an LU factorization given an
%A matrix
%Input
%A - input matrix
%Outputs
%L - Lower triangular matrix
%U - Upper triangular matrix
%p - permutation matrix

m = size(A,1);
n = size(A,2);
status = 1;
L = zeros(m,m);
if m == 1 && n == 1
    L = [1];
    U = A;
    status = 0;
    return;
end
for j = 1:n-1
    for i = j+1:m
        if A(j,j) == 0
            fprintf(2,'Denominator was zero. Had to quit.\n\n');
            U = 0;
            return;
        end
        mult = A(i,j)/A(j,j);
        A(i,:) = A(i,:) - mult*A(j,:);
        L(i,j) = mult;
    end
end

L = L+eye(m);
U = A;
status = 0;
end