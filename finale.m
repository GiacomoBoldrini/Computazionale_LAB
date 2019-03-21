function [u,e,En,x]=finale(a)
% [u,e,En,x]=finale(a) is a program which computes the best Hamiltonian for
% a potential |x^a| and tests them with the quasi-classical approximation 
% of Bohr-Sommerfeld. It also computes Hamiltonians with other methods and 
% compares the results. As outputs it gives the result of the best
% computation of the eigenvectors (u) and eigenvalues (e) using the built
% in function eig, the Bohr-Sommerferld approximation (En) and the lattice
% (x). The program also displays some numerical comparision between the
% eigenfunction in the known case of a=2 (harmonic oscillator).
% If the input value is larger than 113 we exceed realmax in our
% calculations. Since for a->infinity we should get an infinite potential
% well for a>100 another program 'finale1' is called that describe the case
% of an infinite potential well.

if nargin < 1
    a = 2;
end
if  a<=0 
    error('alpha must be a positive scalar')
end
if  ~isscalar(a)
    error('alpha must be a positive scalar')
end
if nargin>1
    error('too many input arguments')
end

if a > 100
    disp('alpha is too big for Matlab to compute. We approximate to potential well')
    finale1;
    return
end

%-------Defining the lattice, wavenumber and potential--------------
N=1024;
Ex= -(N-1)/2:(N-1)/2;
Potential = @(x) abs(x).^a;
V = Potential(Ex);
n = floor(N/2);
nn = floor((N-1)/2);
k = (2*pi/N)*(-n:nn)';

%-------optimal lattice spacing case k.^2/2-------------
maxT=max(k.^2/2);
maxV=max(V);
dx = (maxT/maxV)^(1/(2+a));
x=Ex'*dx;
V = Potential(x);


%-------optimal lattice spacing case 4*sin(k/2).^2------
maxT1 = max(4*sin(k/2).^2);
dx1 = (maxT1/maxV)^(1/(2+a));
x1 = Ex'*dx1;
V1 = Potential(x1);

Z=fft(eye(N));

%-------Kinetic matrix case 4*sin(k/2).^2--------------
k1=fftshift(k); 
T1=ifft(diag(4*(sin(k1/2)).^2)*Z);
T1 = real(T1+T1')/2;
T1(abs(T1)<1e-10)=0;
T1 = T1/dx1/dx1/2;


%-------Kinetic matrix case k.^2/2---------------------
k2=k1/dx;
T=ifft(diag(k2.^2/2)*Z);
T=real(T+T')/2;

%------Tri-Diagonal Laplacian & kinetic matrix----------
eg = ones(N,1);
T2 = full(-0.5*spdiags([eg -2*eg eg],[-1 0 1],N,N) / dx1^2);

%-------eigen values/vectors k.^2/2------------

[u,e]= diagonal(T,V);

%-------eigen values/vectors 4sen(k/2).^2------

[u1,e1]= diagonal(T1,V1);

%-------eigen values 3diag-------------

[~,e3]= diagonal(T2,V1);

%-------Bohr-Sommerfeld quantization----------------

n=0:1:N-1;
En=((pi*(n+0.5))/(2*sqrt(2)*integral(a))).^(2*a/(a+2));

%------------Plotting Eigenvalues--------------------
figure(1)
subplot(1,2,1)
plot(e)
hold on
plot(e1)
plot(e3)
plot(En)
legend('k^2/2','4sen(k/2)^2','Tridiagonal','Sommerfeld')
xlabel('n')
ylabel('e')
title('Energy eigenvalues')
subplot(1,2,2)
plot(abs(e-En'))
legend('|e-E_{n}|')
xlabel('n')
title('k^2 and Bohr-Sommerfeld eigenvalue difference')

%------------Plotting Eigenfunctions--------------------
figure (2)
hold on
plot(x,(u(:,1)))
plot(x,(u(:,2)))
plot(x,(u(:,3)))
legend('u_{1}','u_{2}','u_{3}')
title('eigenfunctions computed with k')
xlabel('x')
ylabel('u')

%------------Eigenfunctions Analysis a=2--------------------
if a==2
figure(3)
a2plots(u1,u,x,dx)
end

%--------plotting probability densities---------
figure(4)
qcprobability(x,u,V,dx,e,N)

end
function [u,e] = diagonal(T,V)
H=T+diag(V);
[u,e]=eig(H);
e=diag(e);
if ~issorted(e)
    [e,I]=sort(e);
    u(:,I);
end
end
function C = integral(a)
syms t;
y = sqrt(1-t.^a);
C = int(y,[0,1]); 
end
