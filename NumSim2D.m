function NumSim2D()
N = 10;
M = 10;
maxit = 50;
T1 = 300 ; T2 = -10 ; 
P1 = 2*T1 ; P2 = 100;
BCP = @(x,y) -x;
cp = 2.08*10^-3; kappa = 0.025; k = 0.01; h_lv = 2.257;
sigma = cp/kappa;
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
    if (max(max(abs(Psol-P)))< 0.01)
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
    figure(2)
    surf(reshape(P,N,M)')
    title('Pressure')
    drawnow
end
[X,Y] = meshgrid(linspace(-1,1,round(N))); 
quiver(ux(X,Y),uy(X,Y));


end