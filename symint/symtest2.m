function symtetst2
%SYMTEST2 Test for symplectic methods.
%  SYMTEST2 computes Poincaré section for
%  
%     dx/dt = -dH/dy, dy/t = -dH/dx 
%
%  with Hamiltonian
%
%     H = sin(pi*x)*sin(pi*y) + epsilon*sin(3*pi*x)*sin(pi*y)*cos(omega*t)
%
%  defined in the unit square using ODE4 (nonadaptive 4th-order RK), ODE45
%  (adaptive Dormand-Prince 5(4) embedded RK pair), GL1 (midpoint rule) and
%  GLS (fully implicit s-stage Gauss-Legendre) methods.

%  Francisco J. Beron-Vera, 22-Mar-2005
%  $Revision: 1.2 $  $Date: 18-May-2005 14:58:31 $

UseDefaultSettings=input('Use default settings? ','s');

if strcmp(UseDefaultSettings(1),'y') || strcmp(UseDefaultSettings(1),'Y')
	epsilon=.25;
	omega=30;
	T=2*pi/omega;
	options=odeset('RelTol',1e-6,'AbsTol',1e-6);
	s=3;
	tol=1e-6;
	maxiter=15;
	Nt=100;
	nt=5;
else
	disp('Parameters of perturbation.')
	epsilon=input('  Strength: ');
	omega=input('  Frequency: ');
	T=2*pi/omega;
	disp('Parameters for ODE45.')
	reltol=input('  Relative tolerance: ');
	abstol=input('  Absolute tolerance: ');
	options=odeset('RelTol',reltol,'AbsTol',abstol);
	disp('Parameters for GLS.')
	s=input('  Stage number: ');
	tol=input('  Tolerance for Newton increments: ');
	maxiter=input('  Maximum Newton iterations: ');
	disp('Integration parameters.')
	Nt=input('  Integration time (multiple of perturbation period): ')/T;
	nt=T/input('  Time step (fraction of perturbation period): ');
end

t=0:T/nt:Nt*T;
ti=T:T:Nt*T;
y=linspace(.51,.99,30); x=repmat(.5,size(y)); z0=[x;y]; z0=z0(:);

I=nt+1:nt:nt*Nt;
Jx=1:2:2*length(x)-1;
Jy=2:2:2*length(x);

clf
subplot(221)
zODE4=ode4(@X,t,z0,epsilon,omega); 
plot(zODE4(I,Jx),zODE4(I,Jy),'.','MarkerS',.1)
axis([0 1 0 1])
title('nonadaptive 4th-order RK')
disp('ODE4 done.')
clear zODE4

subplot(222)
[dummy zODE45]=ode45(@X,[0 ti],z0,options,epsilon,omega);
plot(zODE45(2:end,Jx),zODE45(2:end,Jy),'.','MarkerS',.1)
axis([0 1 0 1])
title('adaptive Dormand-Prince 5(4) embedded RK pair')
disp('ODE45 done.')
clear zODE45

subplot(223)
zGL1=gls(@X,@DX,t,z0,1,ti,tol,maxiter,epsilon,omega); 
plot(zGL1(:,Jx),zGL1(:,Jy),'.','MarkerS',.1)
axis([0 1 0 1])
title('implicit 2nd-order GLRK (midpoint rule)')
disp('GL1 done.')
clear zGL1

subplot(224)
zGLS=gls(@X,@DX,t,z0,s,ti,tol,maxiter,epsilon,omega); 
plot(zGLS(:,Jx),zGLS(:,Jy),'.','MarkerS',.1)
axis([0 1 0 1])
title(['implicit ' num2str(2*s) 'th-order GLRK'])
disp(['GL' num2str(s) ' done.'])
clear zGLS

%--------------------------------------------------------------------------
function out=X(t,z,epsilon,omega)
% Vetor field.
n=length(z);
x=z(1:2:end-1);
y=z(2:2:end);
Hx=pi*cos(pi*x).*sin(pi*y)+...
   epsilon*3*pi*cos(3*pi*x).*sin(pi*y).*cos(omega*t);
Hy=pi*sin(pi*x).*cos(pi*y)+...
   epsilon*  pi*sin(3*pi*x).*cos(pi*y).*cos(omega*t); 
out=[Hy;-Hx];
out=reshape(out,[n/2 2])';
out=out(:); 
%--------------------------------------------------------------------------
function out=DX(t,z,epsilon,omega)
% Jacobian matrix.
x=z(1:2:end-1);
y=z(2:2:end);
Hxx=-pi^2*sin(pi*x).*sin(pi*y)-...
   epsilon*(3*pi)^2*sin(3*pi*x).*sin(pi*y).*cos(omega*t);
Hxy= pi^2*cos(pi*x).*cos(pi*y)+...
   epsilon*3*pi^2  *cos(3*pi*x).*cos(pi*y).*cos(omega*t);
Hyy=-pi^2*sin(pi*x).*sin(pi*y)-...
   epsilon*pi^2    *sin(3*pi*x).*sin(pi*y).*cos(omega*t); 
out=[diag(Hxy) diag(Hyy); -diag(Hxx) -diag(Hxy)];