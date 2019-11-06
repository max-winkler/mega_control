% This script computes the solution of a Poisson equation 
%   -y'' = f
% with randomly chosen Dirichlet conditions on a metric graph. 

clear all;
close all;

addpath('./core/');
addpath('./graphs/');

% Subintervals per edge
ne = 100;

% Create graph object
G = L_graph(5);
m = size(G.Edges,1);    % number of edges
n = size(G.Nodes,1);    % number of nodes

% Randomly select some Dirichlet nodes
nd = floor(n/4);
ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

% Assemble FE matrices and load vector
[Lex,~,F] = assemble(G, ne, 1, 0, 1);

ntil = (ne-1)*m;        % number of interior points overall
nf = ntil+n-nd;         % number of free nodes

FreeNodes = 1:nf;

% Solve system of linear equations
y = zeros(size(Lex,1),1);
y(FreeNodes,1) = Lex(FreeNodes,FreeNodes) \ F(FreeNodes);

% Plot the solution
plot_function_over_graph(G,y,nd);