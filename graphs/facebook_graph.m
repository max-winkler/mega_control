function G = facebook_graph
%FACEBOOK_GRAPH Creates the ego graph of the facebook social network
%   facebook_graph loads the graph data stored in the file facebook.mat and
%   creates an instance of the graph data structure.

% Code snippet is taken from
% https://blogs.mathworks.com/loren/2016/02/03/visualizing-facebook-networks-with-matlab/
load facebook

who
comb = vertcat(edges{:});                           % combine edges
comb = sort(comb, 2);                               % sort edge order
comb = unique(comb,'rows');                         % remove duplicates
comb = comb + 1;                                    % convert to 1-indexing
combG = graph(comb(:,1),comb(:,2));                 % create undirected graph
notConnected = find(degree(combG) == 0);            % find unconnected nodes
G = rmnode(combG, notConnected);

end

