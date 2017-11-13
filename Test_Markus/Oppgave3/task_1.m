function [x,y,U] = task_1(M,N)
%Funksjonen tar inn M for antall aktive x-punkter og N for antall aktive
%y-punkter. Finner så løsningen for laPlace ligningen og plotter den med
%randbetingelser U(x,1) = sin(pi*x)
%% Initialsierer antatt kjente verdier med k=h
fe = @(x,y) (sinh(pi)^-1)*sin(x.*pi)*sinh(y.*pi); %Finn en måte å få denne grafisk
gx1 = @(x) sin(pi*x);
x = linspace(0,1,M+2);
y = linspace(0,1,N+2);
h = x(2)-x(1);
k = y(2)-y(1);
%% Creating A and T
e = ones((M+2)*(N+2),1);
A = spdiags(-2*(k^-2+h^-2)*e,0,M*N,M*N);
e = ones(M,1);
T = spdiags([h^-2*e,-2*(k^-2+h^-2)*e,h^-2*e],-1:1,M ,M);
I = sparse(eye(size(T)));
for i = 1:M:M*N
    A(i:i+M-1,i:i+M-1) = T;
end
for i = 1:M:M*N-M
    A(i:i+M-1,i+M:i+2*M-1) = I*k^-2;
    A(i+M:i+2*M-1,i:i+M-1) = I*k^-2;
end
%% Alternative way
e = ones(M*N,1);
diagonal2 = -2*e*(h^-2+k^-2);
A1 = spdiags([k^-2*e,e*h^-2,diagonal2,e*h^-2,e*k^-2],[-M,-1,0,1,M],M*N,M*N);
indices1 = [M+1:M:M*N]; indices2 = [M:M:M*N-1];
A1(indices1,indices2) = 0; A1(indices2,indices1) = 0;
%% Creating F and U
F = zeros(M*N,1);
F(end-M+1:end) = -(1/k)^2*gx1(x(2:end-1));
Utemp = A\F;
U = zeros(N+2,M+2);
U(end,:) = gx1(x); 
U(2:N+1,2:M+1) = vec2mat(Utemp,M); 
%% Plotter figurer
surf(x,y,U)
xlabel('x')
ylabel('y')
zlabel('U(x,y)')
title('Laplace equation')
figure;
noyaktig = fe(x',y);
surf(x,y,noyaktig');
end


