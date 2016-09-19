function root_table()
%function root_table
%This function calculates all the nth roots of d from d = 2-97 and as n
%ranges from 2-5 via methods of bisection, linear interpolation, secants,
%and Newton
%Inputs
%non
%Outputs
%no formal output but this function does create a table in the console that
%summarizes the results of the entire function

maxitr = 1000;
epsilon = 1e-10;
delta = 1e-10;
loud = 0;

counter = 1;
for i = 2:97
    if isprime(i)
        d(counter) = i;
        counter = counter + 1;
    end
end

n = [2:5];
% f = d^(1/n)
% fp = (1/n)*d^(-1/n);
fprintf('d  n | Bisect.   Lin. Int.   Secants   Newton\n');
fprintf('_____|________________________________________\n');
for i = 1:length(d)
    for j = 1:length(n)
        if isprime(d(i))
            f=@(x) x.^(n(j)) - d(i);
            fp=@(x) (1/n(j))*x^(-1/n(j));
%             fp = (1/n(j))*d(i)^(-1/n(j));
            [xstarBis,fxstarBis, nitrBis, statusBis] = Bisection(f,0,d(i),epsilon,delta,maxitr,loud);
%             if statusBis ~= 0
%                 nitrBis = char(' * ');
%             else
%                 nitrBis = num2str(nitrBis,'%c');
%             end
            [xstarInt, fxstarInt, nitrInt, statusInt] = Interpolation(f,0,d(i),epsilon, maxitr,loud);
%             if statusInt ~= 0
%                 nitrInt = char(' * ');
%             else
%                 nitrInt = num2str(nitrInt,'%c');
%             end
            [xstarSec, fxstarSec, nitrSec, statusSec] = Secant(f,0,d(i),epsilon,maxitr,loud);
%             if statusSec ~= 0
%                 nitrSec = char(' * ');
%             else
%                 nitrSec = num2str(nitrSec,'%c');
%             end
            [xstarNew, fxstarNew, nitrNew, statusNew] = Newton(f,fp,d(i)/2,epsilon,maxitr,loud);
%             if statusNew ~= 0
%                 nitrNew = char(' * ');
%             else
%                 nitrNew = num2str(nitrNew,'%c');
%             end
        end
        if n(j) == 2
            fprintf('%d  %d |  %d     %d     %d      %d\n',d(i),n(j),nitrBis,nitrInt,nitrSec,nitrNew);
        else
            fprintf('    %d |  %d     %d     %d      %d\n',n(j),nitrBis,nitrInt,nitrSec,nitrNew);
        end
    end
    fprintf('________________________________________________\n');
end
fprintf('I couldnt not figure out how to get an asterisk to print in the place of an integer\n');
end