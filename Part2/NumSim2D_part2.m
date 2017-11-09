function NumSim2D_part2()
N = 10;
M = 10;
%% General constants
Qin = 1000; T0 = 200; sigma = 1; gamma = 1; m = 1; h_lv = 1;
%% Vapor constants
Cpv = 1; Kappa_v = 1; Tv = 300; pv = 5; kv = 1; rho_v = 1;
%% Liquid constans
Cpl = 1; Kappa_l = 1; Tl = 300; pl = 5; kl = 1; rho_l = 10;



maxit = 50;
T1 = 300 ; T2 = 200 ; 
BCPv = @(x,y) -x;
BCPl = @(x,y) -x;


Psol = zeros(M*N,1);
ZERO = @(x,y) 0;
for i = 1:maxit
    tic
    [Pv,~,x] = laplace_P_vapor(N, ZERO , BCPv , ZERO , ZERO);             % Solving for P
    [pvx,pvy,BCv] = getGradientALL(Pv,x,N,M);               % Finding the gradient of P
    uvx = @(x,y) -kv*pvx(x,y); uvy = @(x,y) -kv*pvy(x,y);     % Darcy approx
    [Pl,~,x] = laplace_P_vapor(N ,@(x,y) -BCPl(x,y) , ZERO , ZERO ,ZERO);             % Solving for P
    [plx,ply,~] = getGradientALL(Pl,x,N,M);               % Finding the gradient of P
    ulx = @(x,y) -kl*plx(x,y); uly = @(x,y) -kl*ply(x,y);     % Darcy approx
    
    
    
    
    figure(1)
    surf(reshape(T,N,M)')
    title('Temperature')
    drawnow
    figure(2)
    surf(reshape(Pv,N,M)')
    title('Pressure')
    drawnow
end
[X,Y] = meshgrid(linspace(-1,1,round(N))); 
quiver(ux(X,Y),uy(X,Y));


end