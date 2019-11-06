function G = L_graph(n)
%L_GRAPH Creates a graph for a finite difference grid of an L-shaped
%domain.
%   G = L_graph(n) creates a an L-shaped graph where one long edge of the
%   L-shaped domain is divided into n intervals.

A = delsq(numgrid('L',n+3));
G = graph(A,'OmitSelfLoops');

end

