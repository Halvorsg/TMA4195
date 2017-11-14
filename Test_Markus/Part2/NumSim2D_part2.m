function NumSim2D_part2()
addpath ('..\Grids');
addpath ('..\Oppgave1');
N = 20;
M = 20;

%% Vapor constants
Cv_v = 1.46*10^-3; Kappa_v = 2.4; kv = 1; rho_v = 4.85*10^-3; sigma_v = rho_v*Cv_v/Kappa_v;
%% Liquid constans
Cv_l = 4.19*10^-3; Kappa_l = 62; kl = 0.01; rho_l = 1000; sigma_l = rho_l*Cv_l/Kappa_l;
%% General constants
Qin = 1; T0 = 0;  h_lv = 1; 
C_v = 1/(rho_v*kv*h_lv); C_l = 1/(rho_l*kl*h_lv);

maxit = 1;
BCPv =  @(x,y) -x;
BCPl =  @(x,y) -x;
ulx =   @(x,y) -x;
uly =   @(x,y) -x;
T_interface = @(x,y) x-1;



for i = 1:maxit
    tic
    [Pv,~,x] = laplace_P_vapor(N, BCPv);                        % Solving for P
    [pvx,pvy,~] = getGradientALL(Pv,x,N,M);                     % Finding the gradient of P
    uvx = @(x,y) -kv*pvx(x,y); uvy = @(x,y) -kv*pvy(x,y);       % Darcy approx
    
    [Tv,~,x] = temperature_vapor(N,uvx,uvy, @(x,y) Qin, T_interface , T0, sigma_v);
    [~,~,dTv_dy] = getGradientALL(Tv,x,N,M);        % Finding the gradient of P
    Tv = reshape(Tv,M,N)';

    T_interface = @(x,y) interp1(linspace(-1,1,20),Tv(end,:),x);

    [Tl,~,x] = temperature_liquid(N,ulx,uly, @(x,y) Qin, T_interface , T0, sigma_l);
    [~,~,dTl_dy] = getGradientALL(Tl,x,N,M);                    % Finding the gradient of P
       
    BCPl = @(x,y) -(Kappa_l*C_l*dTl_dy(x,y) - Kappa_v*C_l*dTv_dy(x,y));

    [Pl,~,x] = laplace_P_liquid(N ,@(x,y) -BCPl(x,y));          % Solving for P
    [plx,ply,~] = getGradientALL(Pl,x,N,M);                     % Finding the gradient of P
    ulx = @(x,y) -kl*plx(x,y); uly = @(x,y) -kl*ply(x,y);       % Darcy approx
    
    BCPv = @(x,y) (Kappa_l*C_l*dTl_dy(x,y) - Kappa_v*C_v*dTv_dy(x,y));
    
    
    
    
    
    figure(1)
    surf(reshape(Tv,N,M)')
    title('Temperature')
    xlabel('x')
    ylabel('y')
    drawnow
    figure(2)
    surf(reshape(Pv,N,M)')
    title('Pressure')
    drawnow
    figure(3)
    surf(reshape(Tl,N,M)')
    title('Temperature')
    xlabel('x')
    ylabel('y')
    drawnow
    figure(4)
    surf(reshape(Pl,N,M)')
    title('Pressure')
    drawnow
end
figure(3)
[X,Y] = meshgrid(linspace(-1,1,round(N))); 
quiver(uvx(X,Y),uvy(X,Y));


end