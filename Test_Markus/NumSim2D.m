function NumSim2D()
addpath ('..\Grids');
addpath ('..\Oppgave1');
addpath ('C:\Users\halvo\Documents\MATLAB\TMA4195\TMA4195')
N = 30;
M = 30;
maxit = 10;
T1 = 300 ; T2 = -10 ; 
P1 = 2*T1 ; P2 = 100;
BCP = @(x,y) -1*(x)/1.0634;

%https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
%https://en.wikipedia.org/wiki/Enthalpy_of_vaporization
cv = 1.46*10^-3; 
kappa = 2.4; k = 1; h_lv = 2.257; rho_v = 4.85*10^-3;
sigma = rho_v*cv/kappa;
C = kappa/(k*h_lv);

Psol = zeros(M*N,1);
for i = 1:maxit
    tic
    [P,~,x] = laplace_P_Zero(N,BCP,P1,P2);                   % Solving for P
    % Pseduo inverse will set mean(P) = 0
    [px,py,BC] = getGradientALL(P,x,N,M);               % Finding the gradient of P
    ux = @(x,y) -k*px(x,y); uy = @(x,y) -k*py(x,y);     % Darcy approx
    BCT = @(x,y) 1/C*BC(x,y);                           % Interface condition for T
    
    [ T,~,x] = temperature(N,ux,uy,BCT, T1 , T2, sigma);       % Solving for T
    [~,~,BC] = getGradientALL(T,x,N,M);                 % Finding gradient of T
    BCP = @(x,y) C*BC(x,y);                             % Interface condition for P
    disp(max(max(abs(Psol-P))));
    if (max(max(abs(Psol-P)))< 0.0001)
        disp(i)
        disp(max(max(abs(Psol-P))));
        fprintf('Converged')
        break
    end
    Psol = P;
    toc
    figure(1)
    surf(reshape(T,N,M)')
    title('Temperature')
    drawnow
    %figure(2)
    %surf(reshape(P,N,M)')
    %title('Pressure')
    %drawnow
end
n = linspace(-1,1,round(N));
[X,Y] = meshgrid(n); 
h = figure(3);
quiver(X,Y,ux(X,Y),uy(X,Y));
axis tight
y = zeros(1,N);
hold on
plot(n,y-1,'black','linewidth',2)
plot(n,y+1,'black','linewidth',2)
plot(y-1,n,'red','linewidth',2)
plot(y+1,n,'blue','linewidth',2)
title('Velocity field')
saveTightFigure(h,'VelocityField')




end