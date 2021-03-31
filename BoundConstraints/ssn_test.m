% REMARK: This algorithm is implemented as derived in the lecture notes by
% Roland Herzog "Optimierung mit partiellen Differentialgleichungen" from 
% the winter term 2018. See Equation (9.45).

clear all;
close all;

% Include subdirectories
addpath('../core/')
addpath('../graphs/')
addpath('../optimization/')

% Define graph
G = L_graph(10);

m = size(G.Edges,1); % number of edges
n = size(G.Nodes,1);% number of nodes

% Number of intervals per edge
ne = 50;

% Model parameters
mu = 1.; % diffusion
c0 = 1; % potential
src = 0.5;  % source term

nei = ne-1; % number of interior points on an edge
ntil = nei*m; % number of interior points overall

% Define the control vertices
nd = floor(n/3);

ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

nf = ntil+n-nd;

[L,M,f] = assemble(G, ne, mu, c0, src);

yd = M*ones(nf+nd,1);

% Problem parameters
beta = 1e-2;
ub = 1.2;

% Algorithm parameters
param.maxiter = 20;
param.plot = false;
param.verbose = true;

solve_pdas(G,L,M,f,yd,beta,ub,nd,param);
