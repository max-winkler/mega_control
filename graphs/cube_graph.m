function G = cube_graph(n)
%CUBE_GRAPH Creates the finite difference graph of a cube.
%   cube_graph(n) creates the finite difference graph of a cube where each 
%   edge is divided into n intervals.

e=ones(n,1);
D1=spdiags([-e 2*e -e],-1:1,n,n);
I1=speye(n);
A=kron(I1,kron(I1,D1))+kron(D1,kron(I1,I1))+kron(I1,kron(D1,I1));
G = graph(A,'OmitSelfLoops');
end

