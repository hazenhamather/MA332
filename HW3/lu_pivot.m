function [L,U,p] = lu_pivot(A)
%function [L,U,p] = lu_pivot(A)
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
L = zeros(m,m);
p = 1:m;
if m == 1 && n == 1
    L = [1];
    U = A;
    return;
end
if m > n + 1
    limit = n;
elseif n > m+1
    limit = m-1;
elseif m > n
    limit = m - 1;
else 
    limit = n -1;
end
for j = 1:limit
    [colMax,ArowIndex] = max(abs(A(p(j:end),j)));
    if colMax == 0
        fprintf(2,'Denominator was zero. Had to quit.\n\n');
        U = 0;
        return;
    end
    if ArowIndex ~= j
            ArowIndex = ArowIndex + (j-1);
            ptemp = p(j);
            p(j) = p(ArowIndex);
            p(ArowIndex) = ptemp;
    end
    pidx = (p==ArowIndex);
    pidx = find(pidx);
    for i = j+1:m
        mult = A(p(i),j)/A(p(j),j);
        A(p(i),:) = A(p(i),:) - mult*A(p(j),:);
        L(p(i),j) = mult;
    end
end
L = L(p,:);
A = A(p,:);
L = L+eye(m);
U = A;
end
