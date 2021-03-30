clear all;
close all;

addpath('./core/');
addpath('./graphs/');

% Model parameters
eps = 1.e-2;  % diffusion parameter
T = 10;       % Final time

% Subintervals per edge
ne = 20;

% Create graph object
% G = L_graph(5);
G = cube_graph(3);

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

nt=100;
tau = T/nt;
C = spdiags([-ones(nt,1),ones(nt,1)], [0,1], nt-1, nt);

Y = zeros(nall, nt);
Y(FreeNodes,1) = rand(nf,1)-0.5; % Initial value

% Implicit Euler
% for i=1:nt-1
%     
%     fprintf("Time step %d of %d.\n", i, nt);
%     y_new = Y(FreeNodes,i);
%     y_old = Y(FreeNodes,i);
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
%     Y(FreeNodes,i+1) = y_new;
% end

% Newton Method
resid = 1.;
while resid > 1.e-8
    F = Ml(1:nf,1:nf)*Y(1:nf,:)*C' ...
        + tau*eps*K(1:nf,1:nf)*Y(1:nf,2:nt);
    for i=1:nt-1
        F(:,i) = F(:,i) + tau*Ml(1:nf,1:nf)*(Y(1:nf,i+1).^3 - Y(1:nf,i+1));
    end
    F = F(:);        
    
    resid = norm(F);
    fprintf("  newton iteration, resid=%e\n", resid);
    
    if resid > 1.e-8   
        % Solve Newton system        
        %[dy,flag,relres,iter] = minres(@(x)DF(x,K,Ml,C,nall,nf,nt,tau), -F, 1e-8, 10000);
        
        DF = kron(C(:,2:nt), Ml(1:nf,1:nf));
        DF = DF + tau*eps*kron(speye(nt-1,nt-1), K(1:nf,1:nf));
        
        for i=1:nt-1
            DF((i-1)*nf+1:i*nf,(i-1)*nf+1:i*nf) = DF((i-1)*nf+1:i*nf,(i-1)*nf+1:i*nf) ...
                + tau*spdiags(Ml(1:nf,1:nf)*(3*Y(1:nf,i+1).^2-ones(nf,1)),0,nf,nf);
        end
        
        dy = DF\(-F);
        
        % Newton update        
        dY = reshape(dy,nf,nt-1);
        Y(1:nf,2:nt) = Y(1:nf,2:nt) + dY;
    end            
end

% Play video
MOV(nt+1) = struct('cdata',[],'colormap',[]);

for l=1:nt
    plot_function_over_graph(G,Y(:,l));
    title(sprintf("Time step %d of %d", l, nt));
    zlim([min(min(Y)),max(max(Y))]);
    MOV(l) = getframe(gca);
end

movie(MOV)

% Routine for matrix-vector product in Newton system
% function b = DF(y,A,Ml,C,nall,nf,nt,tau)    
%     Y = reshape(y,nall,nt-1);
%     B = Ml(1:nf,1:nf)*Y(1:nf,:)*C(:,2:nt)' + tau*eps*A(1:nf,1:nf)*Y(1:nf,:) + tau*diag(Ml(1:nf,1:nf)).*(3*Y(1:nf,:).^2-1);
%     b = B(:);
% end