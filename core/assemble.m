function [A,M,F] = assemble(G,ne,varargin)
%ASSEMBLE Assembles the finite element stiffness and mass matrix and the
%force vector.
%   [A,M,F] = assemble(G,n,mu,c0,f) assembles   
%       A - the finite element system matrix corresponding to 
%           \int_\Omega (mu u' v'+c0 u v) dx
%       M - the mass matrix corresponding to 
%           \int_\Omega u v dx
%       F - the load vector corresponding to 
%           \int_\Omega f v dx
%   for a given graph G, with ne the number of FE nodes per edge, mu the
%   diffusion parameter, c0 the potential function and f the right-hand
%   side function. Up to now c0 and f have to be a scalar constants. This
%   will be extended in a future version of this code.

m = size(G.Edges,1);    % number of edges

% Set default parameters if not specified
if nargin > 2
    mu = varargin{1};
else
    mu = 1.;
end
if nargin > 3
    c0 = varargin{2};
else
    c0 = 0.;
end

% Create incidence matrices for original nodes
E = -incidence(G);
Ep = 1/2*(E+abs(E));
Em = 1/2*(E-abs(E));

% Create incidence matrices for nodes of extended graph
e1 = zeros(ne,1);e1(1)=1;
en = zeros(ne,1);en(end)=1;
Eh = [];
for j=1:m
    Eh=[Eh,kron(Ep(:,j),e1')+kron(Em(:,j),en')];
end

e = ones(ne,1);
Ee = spdiags([-e e],0:1,ne,ne);Ee(end,:)=[];
Eb = kron(speye(m),Ee);

Et = [Eb;Eh];
he = 1/ne;

% Some auxiliary matrices
We = kron(speye(m),speye(ne)*1/he);
Wm = kron(speye(m),speye(ne)*he);
EWE = abs(Et)*Wm*abs(Et)';

% FE system matrix
L = mu*Et*We*Et';
Mc = c0/6*(abs(Et)*Wm*abs(Et)'+diag(diag(EWE)));
A = L+Mc;

% Mass matrix
M = 1/6*(abs(Et)*Wm*abs(Et)'+diag(diag(EWE)));

% Force vector
if nargin > 4
    F = M*varargin{3}*ones(size(M,1),1);
else
    F = [];
end
