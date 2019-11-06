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

% Lets solve stationary equation on the graph
aux = triu(Lex,1);
%A_control = aux(1:nf,nf+1:end);
A_control = Lex(1:nf,nf+1:end);

A_stiff=Lex; 
A_stiff(:,nf+1:end)=[];
A_stiff(nf+1:end,:)=[];

% I_pen =
beta = 1e-0;

% Create saddle point problem
AA = ...
    [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
    Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
    A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
    ];

% SchurM=Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(Mex(1:nf,1:nf)\Mex(1:nf,nf+1:end)))

b=zeros(nf+nd+nf,1);
b(1:nf+nd,1) = Mex*ones(nf+nd,1);
b(nf+nd+1:end,1) = Mex(1:nf,1:nf)*0*ones(nf,1);

tic
x=AA\b;
toc

%
N=A_control*(diag(diag(Mex(nf+1:end,nf+1:end)+beta*speye(nd)))\A_control');
nn=sum(N,2);
NN=spdiags(sqrt(nn).*sqrt(diag(Mex(1:nf,1:nf))),0,size(N,1),size(N,1));
Khat=(A_stiff+NN);
% Shat=[];
Shat=A_stiff*(diag(diag(Mex(1:nf,1:nf)))\A_stiff)+A_control*((Mex(nf+1:end,nf+1:end)+beta*speye(nd))\A_control');
% Shat=(A_stiff+NN)*(Mex(1:nf,1:nf)\(A_stiff+NN)');
% Shat=A_stiff*(diag(diag(Mex(1:nf,1:nf)))\A_stiff)+A_control*((Mex(nf+1:end,nf+1:end)+beta*speye(nd))\A_control');
% ew=eig(full(blkdiag(Mex(1:nf,1:nf),SchurM)\Mex));

% Solve saddle point system
tic
[x2,flag,relres,iter,resvec] =minres(...
    @(x)matvec_mega(x,Mex,A_stiff,A_control,beta,nd,nf), ...
    b,1e-8,800,...
    @(x)prec_mega(x,Mex,A_stiff,A_control,Khat,Shat,beta,nd,nf));
toc
relres

% hsl_mi20_finalize;

norm(x-x2)/norm(x)

% x=minres(AA,b,100,1e-4);
plot_function_over_graph(G,x(1:nf+nd,1),nd)

return
% Eigenvalue plots
SchurMex=Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(((Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end)))
SchurMapp=Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(diag(diag(Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end)))
PMapp = blkdiag(diag(diag(Mex(1:nf,1:nf))),SchurMapp);
PMex = blkdiag(((Mex(1:nf,1:nf))),SchurMex);
Schurex=[A_stiff A_control]*(PMex\[A_stiff A_control]');
Schurapp=[A_stiff A_control]*(PMapp\[A_stiff A_control]');


N=A_control*(diag(diag(Mex(nf+1:end,nf+1:end)+beta*speye(nd)))\A_control');
nn=sum(N,2);
NN=spdiags(sqrt(nn).*sqrt(diag(Mex(1:nf,1:nf))),0,size(N,1),size(N,1));
Khat=(A_stiff+NN);
Shat=Khat*(diag(diag(Mex(1:nf,1:nf)))\Khat');

Mex2 = [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end);...
    Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd)];

% ew2=eig(full(blkdiag(Mex(1:nf,1:nf),SchurMex)\Mex2));
% plot(sort(real(ew2)),'bo')
% hold on
% ew2=eig(full(blkdiag(Mex(1:nf,1:nf),SchurMapp)\Mex2));
% plot(sort(real(ew2)),'rx')

Pex=blkdiag(PMex,Schurex);
ew=eig(full(Pex\AA));
plot(sort(real(ew)),'rx')
hold on
Papp=blkdiag(PMapp,Schurapp);
ew=eig(full(Papp\AA));
plot(sort(real(ew)),'bo')
Papp=blkdiag(PMapp,Shat);
ew=eig(full(Papp\AA));
plot(sort(real(ew)),'mx')

return
plot_function_over_graph(G,x(1:nf+nd,1),ne,ntil,nf,nd)
