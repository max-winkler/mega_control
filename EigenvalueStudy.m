% This script computes the eigenvalues of the KKT saddle point matrix and
% our preconditioned version.

addpath('./core/')
addpath('./graphs/')

G = L_graph(5);

m = size(G.Edges,1); % number of edges
n = size(G.Nodes,1);% number of nodes

ne = 5; % number of intervals

disp('Assemble FE matrices');

% Randomly select some Dirichlet nodes
nd = floor(n/4);
ind = randperm(n,nd);
ordering=[setdiff(1:n,ind),ind];
G = reordernodes(G,ordering);

[Lex, Mex, ~] = assemble(G, ne, 1, 1);

nei = ne-1; % number of interior points on an edge
ntil = nei*m; % number of interior points overall

nf = ntil+n-nd; % number of non-Dirichlet DOFs

% Assemble KKT matrix
aux = triu(Lex,1);
A_control = aux(1:nf,nf+1:end);
A_stiff=Lex(1:nf,1:nf);

for ii = 1:4 % Number of refinement steps
    
    ne = 2^ii;
    fprintf('\nComputation for ne=%d\n=========================\n', ne)
    
    for jj = 1:4 % Exponent of regularization parameter
        
        beta = 10^(-jj);        
        
        fprintf('   and beta=%e\n', beta)
        
        AA = ...
            [Mex(1:nf,1:nf), Mex(1:nf,nf+1:end), A_stiff;...
            Mex(nf+1:end,1:nf), Mex(nf+1:end,nf+1:end)+beta*speye(nd), A_control';
            A_stiff, A_control, sparse(nf,nf)
            ];        
        
        % Build preconditioners
        SchurMex = Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(((Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end)));
        SchurMapp=Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(diag(diag(Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end)));
        PMapp = blkdiag(diag(diag(Mex(1:nf,1:nf))),SchurMapp);
        PMex = blkdiag(((Mex(1:nf,1:nf))),SchurMex);
        Schurex=[A_stiff A_control]*(PMex\[A_stiff A_control]');
        Schurapp=[A_stiff A_control]*(PMapp\[A_stiff A_control]');
        
        SM=(Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(diag(diag(Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end))));
        N=A_control*(diag(diag(SM))\A_control');
        NN=sqrt(N)*sqrt(((Mex(1:nf,1:nf))));
        Khat=(A_stiff+NN);
        Shat=A_stiff*(Mex(1:nf,1:nf)\A_stiff')+A_control*((Mex(nf+1:end,nf+1:end)+beta*speye(nd)-(Mex(nf+1:end,1:nf)*(diag(diag(Mex(1:nf,1:nf)))\Mex(1:nf,nf+1:end))))\A_control');
        N1 = Mex(1:nf,1:nf);
        N2 = N;
        
        Khat=(A_stiff+NN);
        Shat=Khat*(diag(diag(Mex(1:nf,1:nf)))\Khat');
        Shat2=(A_stiff+N1)*(((Mex(1:nf,1:nf)))\(A_stiff'+N2));
        
        Mex2 = [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end);...
            Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd)];
        
        % ew2=eig(full(blkdiag(Mex(1:nf,1:nf),SchurMex)\Mex2));
        % plot(sort(real(ew2)),'bo')
        % hold on
        % ew2=eig(full(blkdiag(Mex(1:nf,1:nf),SchurMapp)\Mex2));
        % plot(sort(real(ew2)),'rx')
        
        %         Pex=blkdiag(PMex,Schurex);
        %         ew=eig(full(Pex\AA));
        %         plot(sort(real(ew)),'rx')
        %         hold on
        %         Papp=blkdiag(PMapp,Schurapp);
        %         ew=eig(full(Papp\AA));
        %         plot(sort(real(ew)),'bo')
        Papp=blkdiag(PMapp,Shat2);
        
        ew_original = eig(AA);
        ew_preconditioned =eig(full(Papp\AA));
        
        subtitle = sprintf('beta=%e, ne=%d', beta, ne);
        
        figure(1);
        plot(ew_original, 'o')
        title({'Eigenvalues without preconditioning',subtitle})
        
        figure(2)
        plot(ew_preconditioned,'o')
        title({'Eigenvalues with preconditioning',subtitle})
        
        pause;
    end
end
