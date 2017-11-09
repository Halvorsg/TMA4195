function g = getGradient(u,spot,p)
%     [ u,~,spot,p] = temperature(100);
%     plot(p,u(spot));
%% Initializing
    h = diff(p);
    dTdy = zeros(length(spot),1);
%% Using finite differences    
    dTdy(1) = (-3/2*u(spot(1)) + 2*u(spot(2)) - 1/2*u(spot(3)))/h(1);
    dTdy(end) = (3/2*u(spot(end)) + -2*u(spot(end-1)) + 1/2*u(spot(end-2)))/h(2);
    dTdy(2:end-1) = (u(spot(3:end))-u(spot(1:end-2)))/(2*h(1));
%% Creating bounadry function of x and y
    g = @(x,y) spline(p,dTdy,x);  
end