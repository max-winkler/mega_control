function G = star_graph(n)
%STAR_GRAPH Create a simple star-shaped graph
%   G = star_graph(n) returns a graph object with n edges outgoing from a 
%   joint center vertex.

A = sparse(n, n);
A(1,2:end)=1;
A(2:end,1)=1;
G = graph(A);

end

