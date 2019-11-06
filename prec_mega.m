function y=prec_mega(x,Mex,A_stiff,A_control,Khat,Shat,beta,nd,nf)
% PREC_MEGA Schur complement preconditioner for the KKT matrix of the 
% optimal control problem

% block 1,1
y(1:nf,1) = Mex(1:nf,1:nf)\x(1:nf,1);

% block 2,2
y(nf+1:nf+nd,1) = (Mex(nf+1:end,nf+1:end)+beta*speye(nd))\x(nf+1:nf+nd,1);

% block 3,3: Schur-complement block
% y(nf+nd+1:2*nf+nd,1) = Khat\x(nf+nd+1:2*nf+nd,1);
% y(nf+nd+1:2*nf+nd,1) = diag(diag(Mex(1:nf,1:nf)))*y(nf+nd+1:2*nf+nd,1);
% y(nf+nd+1:2*nf+nd,1) = (Khat')\y(nf+nd+1:2*nf+nd,1);
y(nf+nd+1:2*nf+nd,1) = Shat\x(nf+nd+1:2*nf+nd,1);

end

% S=A_stiff*(Mex(1:nf,1:nf)\A_stiff)+A_control*((Mex(nf+1:end,nf+1:end)+beta*speye(nd))\A_control');
% AA = ...
%     [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
%     Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
%     A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
% ];