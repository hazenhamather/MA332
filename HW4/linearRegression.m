function [alpha] = linearRegression(x,y,n)
%function [alpha] = linearRegression(x,y,n)
%This function will compute a vector of alpha values given input data along
%with the length of the alpha vector being n and at most m-1, where m is
%the total number of data points
%Inputs
%x - x-coordinates of the data
%y - y-coordinates of the data
%n - order of the polynomial that will be created with alpha
%Output
%alpha - vector being at most m-1 in length which contains the coefficients
%for the polyval function to approximate the data

x = reshape(x,[],1); %Reshaping x so that no matter what is is one column
y = reshape(y,[],1); %Reshaping y so that no matter what is is one column
[xrow,~] = size(x); %Getting the number of rows in x
[yrow,~] = size(y); %Getting the number of rows in y
alpha = 0; %Default alpha to zero in case of an error

if xrow ~= yrow %Checking that x and y have the same number of rows
    fprintf(2,'x and y need to be the same size\n');
    return;
end

if n >=(xrow -1) %Checking that n is not too big
    fprintf(2,'n must be less than length of x - 1\n');
    return;
end

cols = n+1; %need to have one more column than n
A = zeros(xrow,cols); %initializing our A matrix
for i = 1:xrow %row iterator
    for j = 1:cols %column iterator
        A(i,j) = x(i)^(cols-j); %correcting that alpha needs to be reversed
    end
end
alpha = (A'*A)\(A'*y); %solving for the alpha vector
end