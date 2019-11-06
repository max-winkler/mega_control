function y=matvec_mega(x,Mex,A_stiff,A_control,beta,nd,nf)
% MATVEC_MEGA Computes a matrix-vector product with the KKT matrix of the
% optimal control problem
y(1:nf,1) = Mex(1:nf,1:nf)*x(1:nf,1)+Mex(1:nf,nf+1:end)*x(nf+1:nf+nd,1)+A_stiff*x(nf+nd+1:end,1);
y(nf+1:nf+nd,1) = Mex(nf+1:end,1:nf)*x(1:nf,1)+(Mex(nf+1:end,nf+1:end)+beta*speye(nd))*x(nf+1:nf+nd,1)+A_control'*x(nf+nd+1:end,1);
y(nf+nd+1:2*nf+nd,1) = A_stiff*x(1:nf,1)+(A_control)*x(nf+1:nf+nd,1);

end

% The explicit form of matrix is
% AA = ...
%     [Mex(1:nf,1:nf) Mex(1:nf,nf+1:end) A_stiff;...
%     Mex(nf+1:end,1:nf) Mex(nf+1:end,nf+1:end)+beta*speye(nd) A_control';
%     A_stiff A_control   sparse(size(A_stiff,1),size(A_stiff,1))
% ];