function [ u_sol,u_exact,Neumann_bottom,Neumann_points] = temperature(N)
addpath Grids
addpath Oppgave1
%% Get triangle
[p,tri,edges] = getPlate(N);
%% Functions
u1 = @(x,y) x;
u2 = @(x,y) y;
f = @(x,y) 4*(x^2 - 1)*(y^2 - 1);                   % Right hand side
gtop = @(x,y) 2*(x^2-1);                            % Neumann top
gbottom = @(x,y)  2*(x^2-1);                        % Neumann bottom % Be carefull with the normal vector
exact_solution = @(x,y) (x.^2 - 1).*(y.^2 - 1); % Exact solution
%% Pre-allocating
A = spalloc(length(p),length(p),10*length(p));
M = spalloc(length(p),length(p),10*length(p));
F = zeros(length(p),1);
% Taking care of Neumann boundary
p(abs(p)<100*eps) = 0;      % Setting all values of p approx 0 to 0
edges1 = p(edges(:,1),:);   % Position of first node in edge pair
edges2 = p(edges(:,2),:);   % Position of second node in edge pair
Neumann_edges = zeros(2*N,2); % Pre-allocating table of Neumann edges
Neumann_bottom = unique(edges(p(edges,2) == -1));
Neumann_points = p(Neumann_bottom,1);
%% Neumane edges

% i = 1;
cnt = 0;
for i = 1:length(edges)
% plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'r')
    p1 = edges1(i,:); p2 = edges2(i,:);
    if p1(2) == -1 && p2(2) == -1
            cnt = cnt+1;
            %Vertical <=> x = 0
            Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gbottom);            % Line integral for Neumann boundary
            F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
%             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
            Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
            Neumann_edges2(i,:) = edges(i,:);
            Neumann_tri(i,:) = tri(i,:);
            Neumann_px(i) = p1(1);
            
    elseif p1(2) == 1 && p2(2) == 1
            %Horizontal <=> y = 0
            Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gtop);            % Line integral for Neumann boundary
            F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
%             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
            Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
            %do nothing
    end
end

edges = edges(sum(ismember(edges,Neumann_edges),2) < 2 , : );               % Removing Neumann edges from Dirchlet edges
%% Creating A and finishing F
phi1 = @(x,y,p1,p2,p3) p2(1)*p3(2) - p3(1)*p2(2) - y*(p2(1) - p3(1)) + x*(p2(2) - p3(2));
phi2 = @(x,y,p1,p2,p3) p3(1)*p1(2) - p1(1)*p3(2) + y*(p1(1) - p3(1)) - x*(p1(2) - p3(2));
phi3 = @(x,y,p1,p2,p3) p1(1)*p2(2) - p2(1)*p1(2) - y*(p1(1) - p2(1)) + x*(p1(2) - p2(2));

tic
for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    gradPhi = [p2(2) - p3(2) , p3(2) - p1(2) , p1(2) - p2(2);
               p3(1) - p2(1) , p1(1) - p3(1) , p2(1) - p1(1)]/jacDet;
    Ak = jacDet*(gradPhi'*gradPhi)/2;
% Mass matrix
    phiU1 = @(x,y) [phi1(x,y,p1,p2,p3);phi2(x,y,p1,p2,p3);phi3(x,y,p1,p2,p3)]*u1(x,y);
    phiU2 = @(x,y) [phi1(x,y,p1,p2,p3);phi2(x,y,p1,p2,p3);phi3(x,y,p1,p2,p3)]*u2(x,y);
    
    phi_u1hat = quadrature2D_convection(p1,p2,p3,3,phiU1)/jacDet;
    phi_u2hat = quadrature2D_convection(p1,p2,p3,3,phiU2)/jacDet;
    
    phi_uhat = [phi_u1hat,phi_u2hat];
    Mk = phi_uhat * gradPhi;
%Load elemen vector - Fk
    Fk = quadrature2D(p1,p2,p3,3,f);      
%Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    M(tri(i,:),tri(i,:)) = M(tri(i,:),tri(i,:)) + Mk;
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end
toc
%% Adding and removing boundary
u_sol = zeros(size(p(:,1)));

y = (~ismember(tri,edges)).*tri;
inner_vertices = unique(y);
inner_vertices = inner_vertices(2:end);

%% Lifting and solving
edges_hot = edges(p(edges) == -1);
edges_cold = edges(p(edges) == 1);

gr = zeros(length(p),1);
gr(edges_hot) = 0; % Kelvin
gr(edges_cold) = 0; % Kelvin

G = F- (A+M)*gr;
A = A(inner_vertices,inner_vertices);
M = M(inner_vertices,inner_vertices);
F = G(inner_vertices);
u = (A+M)\F;

u_sol(inner_vertices) = u;
u_sol = u_sol+gr;
u_exact = exact_solution(p(:,1),p(:,2));

%% Plotting
% figure
% trimesh(tri,p(:,1),p(:,2),u_sol)
% s = sprintf('Approximate solution with N = %i', N);
% xlabel('x')
% ylabel('y')
% title(s)
% 
% figure
% trimesh(tri,p(:,1),p(:,2),u_exact)
% title('Exact solution')

figure
trimesh(tri,p(:,1),p(:,2),u_exact-u_sol)
xlabel('x')
ylabel('y')
end