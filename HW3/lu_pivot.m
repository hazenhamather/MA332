function [L,U,p] = lu_pivot(A)
%function [L,U,p] = lu(A)
%This function computes an LU factorization
%input
%A- matrix to be factored
%outputs
%L-lower triangular matrix
%U-upper triangular matrix
%p-permutation matrix

m = size(A,1);
n = size(A,2);
L = zeros(m,m);
p = 1:m;
for j = 1:n-1
%     [colMax,rowIndex] = max(abs(A(j:end,j)));
    idx = find(p == j);
    [colMax,rowIndex] = max(abs(A(p(idx:end),j)));
    rowIndex = rowIndex + (j-1);
    if rowIndex ~= j
        ptemp = p(j);
        p(j) = p(rowIndex);
        p(rowIndex) = ptemp;
        passed = 0;
%         divisor = A(rowIndex,j);
        for i = j:m
            if colMax == 0
                fprintf(2,'Divide by zero error. Exiting');
                U = 0;
                return;
            end
            %         mult = A(i,j)/A(j,j);
            if i ~= p(rowIndex)
                mult = A(i,j)/colMax;
                A(i,:) = -mult*A(p(i),:) + A(p(rowIndex),:);
                if passed == 0
                    L(i+1,j) = mult;
                else
                    L(i,j) = mult;
                end
            else
                passed = 1;
            end
        end
    else
        for i = j+1:m
            if A(j,j) == 0
                fprintf(2,'Divide by zero error. Exiting');
                U = 0;
                return;
            end
            %         mult = A(i,j)/A(j,j);
            mult = A(i,j)/A(j,j);
            A(i,:) = A(i,:) - mult*A(j,:);
            L(i,j) = mult;
        end
        A = A(p,:);
        
    end
end
L = L+eye(m);
A = A(p,:);
U = A;
end
