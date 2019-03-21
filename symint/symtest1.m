function symtest1
%SYMTEST1 Test for symplectic methods.
%  SYMTEST1 Integrates x" + sin(x) = 0 using SEEQ, SEEP, GLS and ODE4 solvers.

%  Francisco J. Beron-Vera, 12-Feb-2005
%  $Revision: 1.3 $  $Date: 20-Jul-2006 01:30:25 $

clf

% Time span.
t=0:.5:200*2*pi;

% Initial conditions.
p0=[-5:-3 .5:.5:2 3:5]';  
q0=repmat(2*pi,size(p0));

% SEEQ
tic, [q p]=seeq(@Xq,@Xp,t,q0,p0); et=toc;
subplot(221)
plot(mod(q,4*pi),p,'.b','MarkerS',1)
axis([pi 3*pi -5 5]), axis square
xlabel('$q$'), ylabel('$p$')
title(['SEEQ (elapsed time = ' num2str(et) ')']) 
disp('SEIQ done.')

% SEEP
tic, [q p]=seep(@Xq,@Xp,t,q0,p0); et=toc;
subplot(222)
plot(mod(q,4*pi),p,'.b','MarkerS',1)
axis([pi 3*pi -5 5]), axis square
xlabel('$q$'), ylabel('$p$')
title(['SEEP (elapsed time = ' num2str(et) ')']) 
disp('SEEP done.')

% GLS
z0=[q0';p0']; z0=z0(:);
tic, z=gls(@X,@DX,t,z0,3,t,1e-3,5); et=toc;
q=z(2:end,1:2:end-1);
p=z(2:end,2:2:end);
subplot(223)
plot(mod(q,4*pi),p,'.b','MarkerS',1)
axis([pi 3*pi -5 5]), axis square
xlabel('$q$'), ylabel('$p$')
title(['GLS (elapsed time = ' num2str(et) ')']) 
disp('GLS done.')

% ODE4
tic, z=ode4(@X,t,z0); et=toc;
q=z(2:end,1:2:end-1);
p=z(2:end,2:2:end);
subplot(224)
plot(mod(q,4*pi),p,'.b','MarkerS',1)
axis([pi 3*pi -5 5]), axis square
xlabel('$q$'), ylabel('$p$')
title(['ODE4 (elapsed time = ' num2str(et) ')']) 
disp('ODE4 done.')

%--------------------------------------------------------------------------
function out=Xq(p)
out=p;
%--------------------------------------------------------------------------
function out=Xp(q)
out=-sin(q); 
%--------------------------------------------------------------------------
function out=X(t,z)
n=length(z);
q=z(1:2:end-1);
p=z(2:2:end);
Hq=sin(q);
Hp=p;
out=[Hp;-Hq];
out=reshape(out,[n/2 2])';
out=out(:);
%--------------------------------------------------------------------------
function out=DX(t,z)
n=length(z);
q=z(1:2:end-1);
Hqq=cos(q);
Hqp=zeros(n/2,1);
Hpp=ones(n/2,1); 
out=[diag(Hqp) diag(Hpp); -diag(Hqq) -diag(Hqp)];
