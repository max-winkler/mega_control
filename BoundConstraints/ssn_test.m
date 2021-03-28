clear all;
close all;
% clc;


% Include subdirectories
addpath('./core/')
addpath('./graphs/')

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
A_stiff = Lex; 
A_stiff(:,nf+1:end)=[];
A_stiff(nf+1:end,:)=[];

% I_pen =
beta = 1e-1;

% Create saddle point problem
AA = ...
    [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
    Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
    A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
    ];

Mu = Mex(nf+1:end,nf+1:end)+beta*speye(nd);
yd = Mex*ones(nf+nd,1);

n = nf+nd;

% Dummies, wird von setup gebraucht.
u_a = -5;
u_b = 5;

alpha = 1e-4;
beta = 1e-2;
tol = 10^-4
maxiter = 20;

iTeraTions =[];

u = zeros(nd,1);

% Evaluating the Integral \int_{\Omega_2(x_1)} p
% A couple of steps of a stationary iteration
steps = 4;
% p_norms = sqrt(I*p);

iter = 0;
display('iter norm_F');
display('===========');
display(sprintf('%2d   %e', iter, norm_F));

resnorm = Inf
while resnorm > tol && iter < maxiter
    iter = iter + 1;
    
	y = A_stiff\(A_control * u);
	p = A_stiff\(Mex(1:nf,1:nf)*(y - yd));

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
    res = u - max(min( -Mu*u*p/alpha, ub), ua);

	J = [...
		Mex(1:nf,1:nf),             sparse(nf,nd), -A_stiff;...
		sparse(nd,nf), spdiags(chi,0,nu,nu)*Mex(nf+1:end,nf+1:end)+beta*speye(nd),  A_control';...
		-A_stiff,             A_control,           sparse(nf,nf)...
		];

	dx = J \ [zeros(ny,1); res; zeros(ny,1)];

    end



