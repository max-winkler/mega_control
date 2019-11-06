clear all;
close all;
% W = gilbert(n); %create example graph
% W = renga(n,0.9,0.3);

%% Create simple star graph from Benzi and Arioli
ne = 10; % number of intervals

G = L_graph(10);
m = size(G.Edges,1);    % number of edges
n = size(G.Nodes,1);    % number of nodes

% Randomly select some Dirichlet nodes
nd = floor(n/4);
ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

% Assemble FE matrices and load vector
[Lex,Mex,F] = assemble(G, ne, 1, 0, 1);

ntil = (ne-1)*m;        % number of interior points overall
nf = ntil+n-nd;         % number of free nodes

FreeNodes = 1:nf;

% create right hand size
b=zeros(size(Lex,1),1);

nt=50;
tau = 1/nt;
x=zeros(size(Lex,1),nt);

x(:,1)=0; % Initial value

% Implicit Euler
for i=1:nt
    x(FreeNodes,i+1)=(Mex(FreeNodes,FreeNodes)+tau*Lex(FreeNodes,FreeNodes)) ...
        \ (Mex(FreeNodes,FreeNodes)*x(FreeNodes,i)+tau*F(FreeNodes));
end

% Play video
M(nt+1) = struct('cdata',[],'colormap',[]);

for l=1:nt+1
    plot_function_over_graph(G,x(:,l));
    zlim([min(min(x)),max(max(x))]);
    M(l) = getframe(gca);
end

movie(M)











