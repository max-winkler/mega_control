function [y,u,p] = solve_pdas(G,L,M,f,yd,beta,ub,nd,param)
%SOLVE_PDAS Solves an optimal control problem with control constraints
%

n_all = size(f,1);
nf = n_all-nd;

% Set default parameters
if nargin < 9
    param = struct;
end
if ~isfield(param, "maxiter") 
    param.maxiter = 20;
end
if ~isfield(param, "plot") 
    param.plot = true;
end
if ~isfield(param, "verbose") 
    param.verbose = true;
end

% Index sets
iF = 1:nf;       % Free nodes
iD = nf+1:nf+nd; % Dirichlet nodes

% Initial guess
u = zeros(nd,1);

% Compute state and adjoint for initial guess
y = L(iF,iF)\(f(iF,1) - L(iF,iD)*u);
p = L(iF,iF)\(yd(iF) - M(iF,iF)*y - M(iF,iD)*u);

% Plot initial iterate
if param.plot
    plot_function_over_graph(G,[y;u],nd);
    pause(1);
end

iter = 0;

nr_changed = 1;
last_chi = zeros(nd,1);

while nr_changed > 0 && iter < param.maxiter
    
    iter = iter + 1;
    
    % Compute new multiplier
    resid = L(iD,iF)*p + M(iD,iF)*y + M(iD,iD)*u - yd(iD) + beta*u;
    mult = -resid;
    
    % Active set
    chi = (mult + beta*(u-ub) > 0);
    
    nr_active = sum(chi);
    nr_changed = sum(abs(chi - last_chi));
    
    last_chi = chi;
    
    I_A = spdiags(chi, 0, nd, nd);
    I_I = speye(nd,nd)-I_A;
    
    % Console output
    if param.verbose
        disp('=============================');
        fprintf("Iteration      : %i\n", iter);
        fprintf("Active nodes   : %i\n", nr_active);
        fprintf("Inactive nodes : %i\n", nd-nr_active);
        fprintf("Active changed : %i\n", nr_changed);
    end
    
    % Assemble Newton system
    % ordering: columns - state, control, adjoint
    %           rows    - adjoint eq., opt. cond, state eq.
    
    % Right-hand side of Newton system
    F = [yd(iF) ; I_I*yd(iD) + ub*I_A*ones(nd,1) ; f(iF)];
    
    % Left-hand side of Newton system
    DF = [M(iF,iF)     , M(iF,iD)                                 , L(iF,iF) ; ...
        I_I*M(iD,iF) , I_I*(M(iD,iD) + beta*speye(nd,nd)) + I_A , I_I*L(iD,iF) ; ...
        L(iF,iF)     , L(iF,iD)                                 , sparse(nf,nf)];
    
    % Solve system of equations
    x = DF \ F;
    
    y = x(iF);
    u = x(iD);
    p = x(nf+nd+iF);
    
    % Plot initial iterate
    if param.plot
        plot_function_over_graph(G,[y;u],nd);
        pause(1);
    end
end
end

