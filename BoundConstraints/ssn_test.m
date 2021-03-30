% REMARK: This algorithm is implemented as derived in the lecture notes by
% Roland Herzog "Optimierung mit partiellen Differentialgleichungen" held 
% the winter term 2018. See Equation (9.45).

clear all;
close all;

% Include subdirectories
addpath('../core/')
addpath('../graphs/')

% Define graph
G = L_graph(5);

m = size(G.Edges,1); % number of edges
n = size(G.Nodes,1);% number of nodes

% Number of intervals per edge
ne = 50;

% Model parameters
mu = 1.; % diffusion
c0 = 1.; % potential
f = 0.;  % source term

nei = ne-1; % number of interior points on an edge
ntil = nei*m; % number of interior points overall

% Define the control vertices
nd = floor(n/3);

ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

nf = ntil+n-nd;

[L,M,f] = assemble(G, ne, mu, c0, f);
M=diag(sum(M,2));

% Lets solve stationary equation on the graph
A_control = L(1:nf,nf+1:end);
A_stiff = L(1:nf, 1:nf);

% Create saddle point problem
% AA = ...
%     [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
%     Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
%     A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
%     ];

yd = M*ones(nf+nd,1);

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
mult = zeros(nd,1);

% Compute state and adjoint for initial guess
y = L(iF,iF)\(f(iF,1) - L(iF,iD)*u);
p = L(iF,iF)\(yd(iF) - M(iF,iF)*y - M(iF,iD)*u);

% Plot initial iterate
plot_function_over_graph(G,[y;u],nd);

iter = 0;

disp('iter norm_F');
disp('===========');
%display(sprintf('%2d   %e', iter, norm_F));

resnorm = Inf;
nr_changed = 1;
last_chi = zeros(nd,1);

while iter < maxiter
    
    iter = iter + 1;    

    % Residuals for all equations
    res_adjoint = L(iF,iF)*p + M(iF,iF)*(y-yd(iF)) + M(iF,iD)*(u-yd(iD));
    res_state   = L(iF,iF)*y + L(iF,iD)*u - f(iF);
    res_control = L(iD,iF)*p + M(iD,iF)*(y-yd(iF)) + M(iD,iD)*(u-yd(iD)) + beta*u + mult;
    res_multip  = -mult+max(0,(u-ub)+mult);
    
    mult = -res_control;
    
	%res = u - max(min(-A_control'*p/alpha, ub), ua);
	%res_norm = sqrt( res' * Mu * res);

    % Active set
	chi = (mult + (u-ub) > 0);
    
	nr_active = sum(chi);
    nr_inactive = nd-nr_active;
	nr_changed = sum(abs(chi - last_chi));

    I_A = spdiags(chi, 0, nd, nd);
    I_I = speye(nd,nd)-I_A;
    
    % Console output
    fprintf("Iteration      : %i\n", iter);
    fprintf("Active nodes   : %i\n", nr_active);
    fprintf("Active changed : %i\n", nr_changed); 
    
    %res = u - max(min( -A_control'*u*p/alpha, ub), ua);

    % Assemble Newton system
    % ordering: columns - state, control, adjoint
    %           rows    - adjoint eq., opt. cond, state eq.
    
    % Right-hand side of Newton system
    F = [yd(iF) ; ub*I_A*ones(nd,1) ; f(iF)];
    
    % Left-hand side of Newton system 
    DF = [M(iF,iF) , M(iF,iD)                             , L(iF,iF) ; ...
          I_I*M(iD,iF) , I_I*(M(iD,iD) + beta*speye(nd,nd)) + I_A, I_I*L(iD,iF) ; ...
          L(iF,iF) , L(iF,iD)                             , zeros(nf,nf)];
    
% 	DF = [...
% 		Mex(1:nf,1:nf),             sparse(nf,nd), -A_stiff;...
% 		sparse(nd,nf), spdiags(chi,0,nu,nu)*Mex(nf+1:end,nf+1:end)+beta*speye(nd),  A_control';...
% 		-A_stiff,             A_control,           sparse(nf,nf)...
% 		];

%	dx = DF \ [zeros(ny,1); res; zeros(ny,1)];
    
    % Solve system of equations
    x = DF \ F;
    
    y = x(iF);
    u = x(nf+iD);
    p = x(nf+nd+iF);

    % Plot initial iterate
    plot_function_over_graph(G,[y;u],nd);
    pause();

end