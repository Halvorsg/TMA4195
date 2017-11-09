function error_evaluation()
addpath C:\Users\halvo\Documents\MATLAB\TMA4195\Project\Oppgave1
addpath C:\Users\halvo\Documents\MATLAB\TMA4195\Project\Oppgave2
addpath C:\Users\halvo\Documents\MATLAB\TMA4195\Project\Oppgave3
addpath C:\Users\halvo\Documents\MATLAB\TMA4195\Project\Grids
values = 2.^(3:6);
% values = [5,10,15,20];
error1 = zeros(size(values));
error2 = zeros(size(values));
tic
for i = 1:length(values)
    [approx,exact] = temperature(values(i));
    N(i) = length(approx);
    error1(i) = max(abs(approx-exact));
    error2(i) = 1/sqrt(length(approx))*norm(approx-exact);
end
toc
loglog(N,error1)
hold on
loglog(N,error2)
loglog(N,N.^-1)

legend('\infty-norm','2-norm','O(h^1)')