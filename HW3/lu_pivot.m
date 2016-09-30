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
oldRowIndex = Inf;
for j = 1:n-1
    %     [colMax,rowIndex] = max(abs(A(j:end,j)));
    %     idx = find(p == j);
    %     [colMax,rowIndex] = max(abs(A(p(idx:end),j)));
    pholder = (p>=j);
    pholder = pholder.*p;
    %         oldRowIndex = p(rowIndex);
    [colMax,rowIndex] = max(abs(A(p(j:end),j)));
    if rowIndex ~= 1
        rowIndex = rowIndex + (j-1);
        if rowIndex ~= p(j)
            ptemp = p(j);
            p(j) = p(rowIndex);
            p(rowIndex) = ptemp;
            pidx = (p==rowIndex);
            pidx = find(pidx);
            passed = 0;
        end
    end
    pidx = (p==rowIndex);
    pidx = find(pidx);
%         divisor = A(rowIndex,j);
        for i = j:m
            if colMax == 0
                fprintf(2,'Divide by zero error. Exiting');
                U = 0;
                return;
            end
            %         mult = A(i,j)/A(j,j);
            if p(i) ~= p(pidx)
                mult = A(p(i),j)/colMax;
                A(p(i),:) = -mult*A(p(pidx),:) + A(p(i),:);
                if passed == 0
                    L(p(i),j) = mult;
%                     L(i+1,j) = mult;
                else
                    L(p(i),j) = mult;
%                     L(i,j) = mult;
                end
            else
                passed = 1;
            end
        end
        
%     else
%         for i = j+1:m
%             if A(j,j) == 0
%                 fprintf(2,'Divide by zero error. Exiting');
%                 U = 0;
%                 return;
%             end
%             %         mult = A(i,j)/A(j,j);
%             mult = A(i,j)/A(j,j);
%             A(i,:) = A(i,:) - mult*A(j,:);
%             L(i,j) = mult;
%         end

%         A = A(p,:);
        
%     end
    oldRowIndex = (p==rowIndex);
    oldRowIndex = find(oldRowIndex);
end
L = L(p,:);
L = L+eye(m);
A = A(p,:);
U = A;
end
