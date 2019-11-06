% This script determines numerically the convergence rates of the
% discretization of a Dirichlet control problem on a metric graph.

clear all;
close all;

addpath('./core/');
addpath('./graphs/');

G = L_graph(15);

m = size(G.Edges,1); % number of edges
n = size(G.Nodes,1);% number of nodes

% Number of mesh refinement steps
NrComp = 10;

% Number of Dirichlet vertices
nd = 20;

% Regularization parameter
beta = 1.e-2;

% Randomly select some Dirichlet nodes
ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

% Allocate vectors containing the solution for all refinement steps
u = zeros(nd,NrComp);
y = {};
err_state_l2 = zeros(1,NrComp-1);
err_state_h1 = zeros(1,NrComp-1);
err_control = zeros(1,NrComp-1);
 
for ne=2.^(1:NrComp)
    
    ItComp = log2(ne);
    
    fprintf("Computing refinement level %d/%d\n", ItComp, NrComp)
    
    nei = ne-1; % number of interior points on an edge
    ntil = nei*m; % number of interior points overall    
        
    nf = ntil+n-nd;             
            
    [Lex,Mex,F] = assemble(G, ne, 1, 1, 1.5);
            
    % Create saddle point problem
    % Ordering of equations: adjoint eq., opt. condition, state eq
    % Ordering of variables: state, control, adjoint state
    AA = [...
        Mex(1:nf,1:nf), Mex(1:nf,nf+1:end), Lex(1:nf,1:nf)';...
        Mex(nf+1:end,1:nf), Mex(nf+1:end,nf+1:end)+beta*speye(nd), Lex(1:nf,nf+1:end)';
        Lex(1:nf,1:nf), Lex(1:nf,nf+1:end), sparse(nf,nf)
        ];
    
    b=zeros(nf+nd+nf,1);
    % Incorporate rhs yd=1
    b(1:nf+nd,1) = Mex*ones(nf+nd,1);
    % Incorporate rhs f
    b(nf+nd+1:end,1) = F(1:nf,1);
    
    % Solve equation system with the best solver we have
    x=AA\b;

    % Extract control and state
    u(:,ItComp) = x(nf+1:nf+nd);
    y{ItComp} = x(1:nf+nd);
    
    % Plot the solution
    plot_function_over_graph(G,y{ItComp},nd);
    
    %disp('Press button to continue ...')
    %pause
    pause(0.5)
end

% Compute error of the state variable (by comparison with solution on the
% finest grid)
for k=1:NrComp-1
    diff = u(:,k)-u(:,NrComp);
    err_control(k) = norm(diff);
    
    % Prolongate the state variable
    coarse = y{k};
    for i=k:NrComp-1
        fine = prolongate(G, coarse, 2^i);    
        coarse = fine;
    end
        
    % Compute error norm
    diff = fine-y{NrComp};    
    err_state_l2(k) = sqrt(diff'*Mex*diff);
    err_state_h1(k) = sqrt(diff'*Lex*diff);
end

H = 2.^(-[1:NrComp-1]);

% Compute convergence rates
eoc_control = log(err_control(2:end)./err_control(1:end-1))/log(1/2);
eoc_state_l2 = log(err_state_l2(2:end)./err_state_l2(1:end-1))/log(1/2);
eoc_state_h1 = log(err_state_h1(2:end)./err_state_h1(1:end-1))/log(1/2);

% Table for errors and convergence rates
T = table(H(2:end)', err_control(2:end)', eoc_control', ...
    err_state_h1(2:end)', eoc_state_h1', ...
    err_state_l2(2:end)', eoc_state_l2');
T.Properties.VariableNames = {'h', '|u-uh|', 'eoc(u)', ...
    '|y-yh|_H1', 'eoc(y,H1)', ...
    '|y-yh|_L2', 'eoc(y,L2)'};
disp(T)

% Convergence plot
figure(2);
H = 2.^(-[1:NrComp-1]);
loglog(H,err_control,'*-',H,err_state_h1,'o-',H,err_state_l2,'d-');
legend('|u-uh|_2','|y-yh|_H1','|y-yh|_L2');
xlabel('h');
ylabel('err');

