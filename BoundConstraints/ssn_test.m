clear all;
close all;
% clc;


% This is the implementation of a distributed control problem. Here
% M1=M2=N=M (mass matrix). The variable 4_20 refers to an example with 4
% refinements in space and 20 time-steps in time. The other variables work
% in an analogous manner.
% The variable r1 is the matrix from of the observed state in the form
% My_obs.
% cd ../
% cd data_stokes
load M_Stokes.dat
load A_Stokes.dat
load B_Stokes.dat
load Mp_Stokes.dat
load Kp_Stokes.dat
load Mybar_Stokes.dat
load rhs_Stokes.dat
M=spconvert(M_Stokes);
A=spconvert(A_Stokes);
B=spconvert(B_Stokes);
Mp=spconvert(Mp_Stokes);
Kp=spconvert(Kp_Stokes);
r1 = spconvert(Mybar_Stokes)';
r3 = spconvert(rhs_Stokes)'; 
cd ..
% cd('/afs/mpi-magdeburg.mpg.de/home/stollm/from_scratch/paper/In Preparation/sparsity/matlab')
nv = size(M,1); % Matrix size, dimension of spatial discretization
np = size(Kp,1);

Kp(:,1) = 0;
Kp(1,:) = 0;
Kp(1,1) = 1;
n = nv+np;
% y_d_handle = @(x,y)(sin(2*pi*x).*sin(2*pi*y).*exp(2*y)/6);
% y_d_handle = @(x,y)(sin(2*pi*x).*sin(2*pi*y).*cos(2*x*y)/6);
% Dummies, wird von setup gebraucht.
u_a = -5;
u_b = 5;

alpha = 1e-4;
beta = 1e-2;

% Um auf Georgs Abbruchbedingung zu kommen
epsilon = 10^-8 / alpha * beta
maxiter = 20;

iTeraTions =[];

% Constraint on 2-norm of u on stripes:
% u_b = 18;
u_b = +1e10;

epsilon = 1e-6;

u = zeros(nv,1);

y_dr = r1(:,1);
bnd  = r3(:,1);


% Evaluating the Integral \int_{\Omega_2(x_1)} p
% A couple of steps of a stationary iteration
b = [M*u;zeros(np,1)]-bnd;
t=zeros(nv+np,1);
steps = 4;
for i=1:steps
    t=t+blkdiag(A,Mp)\(b-[A B';B sparse(np,np)]*t);
end
yv = t(1:nv,1);
b=([y_dr;zeros(np,1)] - blkdiag(M,sparse(np,np))*t);
t=zeros(nv+np,1);
for i=1:steps
    t=t+blkdiag(A,Mp)\(b-[A B';B sparse(np,np)]*t);
end
pv = t(1:nv,1);
% I = kron(ones(breite2,breite2), speye(breite1,breite1))*h * spdiags(p,[0],breite1*breite2,breite1*breite2); % Mass matrix here for FEM?
% p_norms = sqrt(I*p);
F = u - 1/alpha .* max(0,pv-beta)-1/alpha.*min(0,pv+beta);
norm_F = sqrt(F' * M * F);

iter = 0;
display('iter norm_F');
display('===========');
display(sprintf('%2d   %e', iter, norm_F));

h=sqrt(min(eig(M)));
control = hsl_mi20_control;
control.coarse_solver = 1;
control.v_iterations = 2;
inform = hsl_mi20_setup(Kp,control);


% control2 = hsl_mi20_control2;
% control2.coarse_solver = 2;
% control2.v_iterations = 2;
% control2.pre_smoothing = 5;
inform2 = hsl_mi20_setup2(A+1/sqrt(alpha)*M,control);    


while norm_F > epsilon && iter < maxiter
    iter = iter + 1;

	chi1 = (pv+beta <=0);
	chi2 = (pv-beta >=0);
    % Creating the G matrix
    G= spdiags(chi1+chi2,[0],nv,nv);
    % Creating the right hand side
%     rhs_better = [-[y_dr;zeros(np,1)];-M*F+M*u-1/alpha*G*M*pv;bnd];
    rhs_better = [zeros(nv+np,1);-M*F;zeros(nv+np,1)];
    % Call GMRES
    [x_better,flag,relres,iter_gmres,resvec] = gmres(@(x)Amult_stokes(x,M,A,B,G,nv,np,alpha),rhs_better,[],1e-6,50,@(x)Aprec_stokes(x,M,A,B,Mp,Kp,G,nv,np,alpha,h));
%     x_better=AA\rhs_better;   
    iTeraTions(iter) = iter_gmres(2)

    u  = u +x_better(n+1:n+nv);
    pv = pv+x_better(n+nv+1:n+2*nv);
    yv = yv+x_better(1:nv);
    F = u - 1/alpha .* max(0,pv-beta)-1/alpha.*min(0,pv+beta);
    norm_F = sqrt(F' * M * F);
	display(sprintf('%2d   %e', iter, norm_F));
end
hsl_mi20_finalize2
% load('ind_5_100.dat')
% surf(vecplotdealii(u,ind_5_100,5));
% surf(vecplotdealiiCD(u+bndsol(:,1),index_CD,6));
% colormap hot
% Plot u:
% mesh(reshape(u, sqrt(n), sqrt(n)))
% % Norm of u on stripes:
% sqrt(h*sum(reshape(u, breite1, breite2).^2,2));



