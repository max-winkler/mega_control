clear all;
close all;

% Include subdirectories
addpath('../core/')
addpath('../graphs/')

% Define graph
G = L_graph(20);

m = size(G.Edges,1); % number of edges
n = size(G.Nodes,1);% number of nodes

% Number of intervals per edge
ne = 50;

% Model parameters
mu = 1.; % diffusion
c0 = 1.; % potential
f = 1.;  % source term

nei = ne-1; % number of interior points on an edge
ntil = nei*m; % number of interior points overall

% Define the control vertices
nd = floor(n/3);
ind = randperm(n,nd);

nf = ntil+n-nd;

NN=numnodes(G);
ordering=[setdiff(1:NN,ind),ind];
G = reordernodes(G,ordering);

[Lex,Mex,F] = assemble(G, ne, mu, c0, f);
Mex=diag(sum(Mex,2));

% Lets solve stationary equation on the graph
A_control = Lex(1:nf,nf+1:end);
A_stiff = Lex(1:nf, 1:nf);

%A_stiff(:,nf+1:end)=[];
%A_stiff(nf+1:end,:)=[];

% I_pen =
beta = 1e-1;

% Create saddle point problem
% AA = ...
%     [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
%     Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
%     A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
%     ];

Mu = Mex(nf+1:end,nf+1:end)+beta*speye(nd);
yd = Mex*ones(nf+nd,1);

n = nf+nd;

% Index sets
iF = 1:nf;       % Free nodes
iD = nf+1:nf+nd; % Dirichlet nodes

% Control bounds
ua = -5;
ub = 5;

% Problem parameters
alpha = 1e-4;
beta = 1e-2;

% Algorithm parameters
tol = 10^-4;
maxiter = 20;

% Initial guess
u = zeros(nd,1);
mu = zeros(nd,1);

% Compute state and adjoint for initial guess
y = Lex(iF,iF)\(F(iF,1) - Lex(iF,iD)*u);
p = Lex(iF,iF)\(Mex(iF,iF)*(y - yd(iF)) + Mex(iF,iD)*(u - yd(iD)));

iter = 0;
display('iter norm_F');
display('===========');
%display(sprintf('%2d   %e', iter, norm_F));

resnorm = Inf;

while resnorm > tol && iter < maxiter
    iter = iter + 1;    

	res = u - max(min( -A_control'*p/alpha, ub), ua);
	res_norm = sqrt( res' * Mu * res);

	chi = (ua < -A_control'*p/alpha) & (-A_control'*p/alpha < ub);
	nr_active = sum(chi);
	if j > 0
		nr_changed = sum(abs( chi - last_chi ));
	else
		nr_changed = 0;
	end
	fprintf('% 3d  %e (%d, %d)', j, res_norm, nr_active, nr_changed);
    res = u - max(min( -A_control'*u*p/alpha, ub), ua);

	J = [...
		Mex(1:nf,1:nf),             sparse(nf,nd), -A_stiff;...
		sparse(nd,nf), spdiags(chi,0,nu,nu)*Mex(nf+1:end,nf+1:end)+beta*speye(nd),  A_control';...
		-A_stiff,             A_control,           sparse(nf,nf)...
		];

	dx = J \ [zeros(ny,1); res; zeros(ny,1)];

end



