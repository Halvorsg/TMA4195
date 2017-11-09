function [ux,uy,BC] = getGradientALL(u,p,N,M)
%     [ u,~,spot,p] = temperature(100);
%     plot(p,u(spot));
%     f = @(x,y) -5*x.^3+-y.^3.*x.^2+y.^2+1;
%     fx = @(x,y) -15*x.^2-2*y.^3.*x;
%     fy = @(x,y) -3*x.^2.*y.^2 + 2*y;
%     u = reshape(u,M,N)
%     n = linspace(-1,1,N);
%     [X,Y] = meshgrid(n);
%     u = f(X,Y);
%     surf(n,n,u)
%     figure
%     hold on
    u = reshape(u,N,M)';
    [ux,uy] = gradient(u,2/N,2/M);
    x = linspace(-1,1,N);
    y = linspace(-1,1,M);
    [X,Y] = meshgrid(x,y);
    ux = @(x,y) interp2(X,Y,ux,x,y);
    uy = @(x,y) interp2(X,Y,uy,x,y);
    BC = @(x,y) uy(x,-1);
% %% Initializing
%     h = diff(p);
%     dTdy = zeros(M,N);
%     dTdx = zeros(M,N);
% %% Using finite differences    
%     dTdy(1,:) = (-3/2*u(1,:) + 2*u(2,:) - 1/2*u(3,:))/h(1);
%     dTdy(end,:) = (3/2*u(end,:) + -2*u(end-1,:) + 1/2*u(end-2,:))/h(2);
%     dTdy(2:end-1,:) = (u(3:end,:)-u(1:end-2,:))/(2*h(1));
%     
%     dTdx(:,1) = (-3/2*u(:,1) + 2*u(:,2) - 1/2*u(:,3))/h(1);
%     dTdx(:,end) = (3/2*u(:,end) + -2*u(:,end-1) + 1/2*u(:,end-2))/h(2);
%     dTdx(:,2:end-1) = (u(:,3:end)-u(:,1:end-2))./(2*h(1));
%     
%     BC = @(x,y) spline(p,dTdy(1,:),x);
% %% Creating bounadry function of x and y
%     [X,Y] = meshgrid(p);
%     uy = @(x,y) interp2(X,Y,flip(dTdy),x,y);
%     ux = @(x,y) interp2(X,Y,flip(dTdx),x,y);  
end