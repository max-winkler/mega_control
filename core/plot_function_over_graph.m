function plot_function_over_graph(G,u,nd)
% PLOT_FUNCTION_OVER_GRAPH Plots a given function over a given graph
% network.
%   plot_function_over_graph(G,u,nd) plots the function u over
%   the graph G. The number nd is the number of Dirichlet nodes.
%
%   plot_function_over_graph(G) plots the graph only.

if nargin < 3
    nd = 0;    
end

Edges=table2array(G.Edges(:,1)); % Get Edges from graph
NrEdges = size(Edges,1); % How many egdes are there?

% Create the force layout of the graph
G = matlab.internal.graph.MLGraph(Edges(:,1), Edges(:,2), max(Edges(:)));
[x,y] = G.subspaceLayout(NrEdges,2);
[x,y] = G.forceLayout(x, y, 100, [], 'off');

NrNodes = size(x,1);

% Compute number of intervals per edge
if nargin > 1
    ne = 1+(size(u,1)-NrNodes)/NrEdges;
    if round(ne) ~= ne
        disp(strcat('WARNING: The size of the input vector u is not ', ...
            'reasonable. Please check the dimension.'));
        return;
    end
else
    ne = 1;
end

% Number of interior vertices
ntil = (ne-1)*NrEdges;

% Number of non-Dirichlet nodes
nf = ntil+NrNodes-nd;

% Plotting
for k=1:NrEdges
    xs = linspace(x(Edges(k,1)), x(Edges(k,2)), ne+1);
    ys = linspace(y(Edges(k,1)), y(Edges(k,2)), ne+1);    
    
    if nargin > 1
        zs = [u(ntil+Edges(k,1),1) u((k-1)*(ne-1)+1:k*(ne-1),1)' u(ntil+Edges(k,2),1) ];
    
        % Plot function
        plot3(xs,ys,zs,'b-','LineWidth',3);
        hold on;
    end
    
    % Plot graph
    plot3(xs,ys,zeros(size(xs)),'r-','LineWidth',2);
    hold on;

    %     plot3([xs(1) xs(1)], [ys(1) ys(1)], [0 zs(1)], 'g-')
    %     plot3([xs(end) xs(end)], [ys(end) ys(end)], [0 zs(end)], 'g-')   
    %     pause
end

% Draw values in Dirichlet nodes
if nargin > 2
    plot3(x(end-nd+1:end),y(end-nd+1:end),u(nf+1:end),'go','LineWidth',2);

    % Draw connection between node and function value
    for k=1:NrNodes
        plot3([x(k) x(k)],[y(k) y(k)], [0,u(ntil+k)],'b:','LineWidth',1)
    end
end

hold off;

% Turn off axes. 
% set(gca,'visible','off')

    