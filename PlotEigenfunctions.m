% In this test script an the eigenvalue distribution and the eigenvectors
% of a simple PDE on a graph network will be illustrated.

clear all;
close all;

addpath('./core/');
addpath('./graphs/');

m = 10; % number of edges
n = m+1; % number of nodes
ne = 20; % number of intervals

G = star_graph(n);

[Lex,Mex,~] = assemble(G, ne, 1, 0);

[EV,EW]=eigs(Lex,Mex,size(Lex,1));
[EW,ind]=sort(diag(EW),'ascend');

figure(1);
plot(EW,'*');
title('Eigenvalues')

for k=1:9
    u=EV(:,k);

    nei = ne-1;

    figure(2);
    subplot(3,3,k);
    plot_function_over_graph(G, u);
    title(sprintf('Eigenfunction %d',k));
end









