function [z] = cubicSpline(x,y,v)
%function [z] = cubicSpline(x,y,v)
%This function creates a cubic spline in order to approximate the data
%given
%Inputs
%x - x-coordinates of the data
%y - y-coordinates of the data
%v - an optional vector 
%Output
%z - spline data

z = []; %Initializing z matrix for return, dependent upon v
plot(x,y,'o'); %plotting the given data points
hold on
x = reshape(x,[],1); %making sure the input will be a column vector
y = reshape(y,[],1); %making sure the input will be a column vector
nx = length(x); %getting the length of input x for later use
deltax = zeros((nx-1),1); %initializing the delta-x vector
for i = 1:length(deltax) %iterating through the delta-x vector
    deltax(i) = x(i+1) - x(i); %Sstting the delta-x step
end
M = zeros(nx-2); %initializing system for s variables
u = zeros(nx-2,1); %initializing system for s variables
M(1,1) = 2*(deltax(2)+deltax(1)); 
M(1,2) = deltax(2); 
M(nx-2,nx-2)=2*(deltax(nx-2)+deltax(nx-1)); 
M(nx-2,nx-3)=deltax(nx-2); 
u(1)=6*(((y(3)-y(2))/(deltax(2)))-((y(2)-y(1))/deltax(1))); 
u(nx-2)=6*(((y(nx)-y(nx-1))/(deltax(nx-1)))-((y(nx-1)-y(nx-2))/deltax(nx-2)));

for j=2:(nx-2) %populating matrices to solve for s vector
    M(j,j)=2*(deltax(j+1)+deltax(j));
    M(j,j+1)=deltax(j+1);
    M(j,j-1)=deltax(j);
    u(j)=6*(((y(j+2)-y(j+1))/(deltax(j+1)))-((y(j+1)-y(j))/deltax(j)));
end
s=M\u;
s=[0;s;0];
a=zeros(nx,1);
b=zeros(nx,1);
c=zeros(nx,1);
d=zeros(nx,1);
fullSpline = {};
for i= 1:(nx-1)
    a(i)=y(i);
    b(i)=(((y(i+1)-y(i))/deltax(i))-(((s(i+1)+2*s(i))*deltax(i))/6));
    c(i)=(s(i)/2);
    d(i)=((s(i+1)-s(i))/(6*deltax(i)));
    w = linspace(x(i), x(i)+deltax(i),1000);
    spline=a(i)+b(i).*(w-x(i))+c(i).*(w-x(i)).^2+d(i).*(w-x(i)).^3;
    plot(w,spline);
    fullSpline{i} = @(k)a(i)+b(i)*(k-x(i))+c(i)*(k-x(i)).^2+d(i)*(k-x(i)).^3;
end
fullSpline = fullSpline';
if exist('v')
    z = zeros(1,length(v));
    for i = 1:length(v)
        splineNumber = 0;
        for j = 1:nx
            if v(i) <= x(j)
                z(i) = fullSpline{splineNumber}(v(i));
                break;
            else
                splineNumber = splineNumber + 1;
            end
        end
    end
end
end