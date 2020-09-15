clear all;
close all;

addpath('./core/');
addpath('./graphs/');

% Model parameters
eps = 1.e-2;  % diffusion parameter
T = 20;       % Final time

% Subintervals per edge
ne = 50;

% Create graph object
G = L_graph(5);
m = size(G.Edges,1);    % number of edges
n = size(G.Nodes,1);    % number of nodes

% Randomly select some Dirichlet nodes
nd = floor(n/4);
nd = 0;
ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

% Assemble FE matrices and load vector
[K,M,~] = assemble(G, ne, 1, 0, 0);

nall = size(K,1);       % number of all nodes
ntil = (ne-1)*m;        % number of interior points overall
nf = ntil+n-nd;         % number of free nodes

% Mass lumping
Ml_diag = M*ones(nall,1);
Ml = spdiags(Ml_diag,0,nall,nall);

FreeNodes = 1:nf;

nt=50;
tau = T/nt;

y = zeros(nall, nt);
y(FreeNodes,1) = rand(nf,1)-0.5; % Initial value

% Implicit Euler
% for i=1:nt-1
%     
%     fprintf("Time step %d of %d.\n", i, nt);
%     y_new = y(FreeNodes,i);
%     y_old = y(FreeNodes,i);
%     
%     % Newton method
%     resid = 1.;
%     while resid > 1.e-8
%         
%         % Right-hand side of Newton system
%         Fn = Ml(FreeNodes,FreeNodes)*(y_new - y_old) ...
%             + tau*eps*K(FreeNodes,FreeNodes)*y_new ...
%             + tau*Ml_diag(FreeNodes).*(y_new.^3-y_new);
%                 
%         resid = norm(Fn);
%         
%         fprintf("  resid for Newton method: %e\n", resid);
%         
%         if resid > 1.e-8
%             % System matrix of Newton system
%             DFn = Ml(FreeNodes,FreeNodes) ...
%                 + tau*eps*K(FreeNodes,FreeNodes) ...
%                 + tau*spdiags(Ml_diag(FreeNodes).*(3*y_new.^2-1),0,nf,nf);
%         
%             % Compute Newton step
%             dy = DFn \ (-Fn);
%         
%             y_new = y_new + dy;
%         end
%     end
%         
%     y(FreeNodes,i+1) = y_new;
% end

% Routine for matrix-vector product in Newton system
function b = DF(y,A,Ml,nall,nf,nt)
    Y = reshape(y,nall,nt);
    B = Ml(1:nf)*Y(1:nf,:)*C' + tau*eps*A(1:nf,:)*Yh(1:nf,:) + tau*diag(Ml(1:nf)).*(3*Y(1:nf,:).^2-1);
    b = B(:);
end

% Newton Method
resid = 1.
while resid > 1.e-8        
    F = [kron([ones(nt,1)],tau*M11*ones(nf,1));kron([ones(nt,1)],tau*M22*ones(nd,1));zeros(nt*nf,1)];
end

% Play video
MOV(nt+1) = struct('cdata',[],'colormap',[]);

for l=1:nt
    plot_function_over_graph(G,y(:,l));
    title(sprintf("Time step %d of %d", l, nt));
    zlim([min(min(y)),max(max(y))]);
    MOV(l) = getframe(gca);
end

movie(MOV)