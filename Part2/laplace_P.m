function [ u_sol,Neumann_bottom,Neumann_points] = laplace_P(N,gbottom, P1,P2)
addpath Grids
addpath Oppgave1
%% Get triangle
[p,tri,edges] = getPlate(N);
% TR = triangulation(tri, p);
% triplot(TR)
% hold on
%% Functions
f = @(x,y)      0;                          % Right hand side
gtop = @(x,y)   0;                          % Neumann top
% gbottom = @(x,y) -gbottom(x,y);
fcn = @plus;
%% Pre-allocating
A = spalloc(length(p),length(p),10*length(p));
F = zeros(length(p),1);
phi = [eye(2),[-1;-1]];
%% Taking care of Neumann boundary
p(abs(p)<100*eps) = 0;      % Setting all values of p approx 0 to 0
edges1 = p(edges(:,1),:);   % Position of first node in edge pair
edges2 = p(edges(:,2),:);   % Position of second node in edge pair
Neumann_edges = zeros(2*N,2); % Pre-allocating table of Neumann edges

Neumann_bottom = unique(edges(p(edges,2) == -1));
Neumann_points = p(Neumann_bottom,1);

%% Neumane edges

% i = 1;
for i = 1:length(edges)
% plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'r')
    p1 = edges1(i,:); p2 = edges2(i,:);
    if p1(2) == -1 && p2(2) == -1
            %Vertical <=> x = 0
            Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gbottom);            % Line integral for Neumann boundary
            F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
%             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
            Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
            
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
parfor i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/2;
%Load elemen vector - Fk
       
%Adding into Stiffnessmatrix and load vector
    [rowMap,colMap] = get_node_mapping_matrix(tri(i,:));
    A = fcn(A,sparse(rowMap,colMap,Ak,length(p),length(p)));
    
end

for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
    Fk = quadrature2D(p1,p2,p3,3,f);   
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end
%% Adding and removing boundary
u_sol = zeros(size(p(:,1)));

y = (~ismember(tri,edges)).*tri;
inner_vertices = unique(y);
inner_vertices = inner_vertices(2:end);

%% Lifting and solving
edges_hot = edges(p(edges) == -1);
edges_cold = edges(p(edges) == 1);

gr = zeros(length(p),1);
gr(edges_hot) = P1; % Kelvin
gr(edges_cold) = P2; % Kelvin

G = F - A*gr;
A = A(inner_vertices,inner_vertices);
F = G(inner_vertices);
u = A\F;

u_sol(inner_vertices) = u;
u_sol = u_sol+gr;

%% Plotting
% figure
% trimesh(tri,p(:,1),p(:,2),u_sol)
% s = sprintf('Approximate solution with N = %i', N);
% title(s)
end
    
    
    
    
    