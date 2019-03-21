function finale1
%finale1 computes eigenvalues and eigenvectors for a symmetric (-1,1)
%infinite potrential well. It is called by 'finale' when the input
%parameter is too large for Matlab to compute.

%--------Particle in a box (-1,1)----------

L=1;
N=1024;
dx = 2*L/N;
x  = -L*(1-1/N) : dx : L*(1-1/N);
k = (pi/dx/N)*(1:N)';

%---------kinetic energy for particle in  a box----
T = sinft(diag(k.^2/2)*sinft(eye(N)));
T = real(T+T')/2;
    
%-------Eigenvalues  particle in  a box----
[u,e]=eig(T);
e=diag(e);
if ~issorted(e)
    [e,I]=sort(e);
    u(:,I);
end

%-------Exact eigenvalues---------
En =((1:N).^2)*pi^2/(2*(dx*N)^2);

%-------plotting eigenvalues---------
figure(1)
subplot(1,2,1)
plot(e)
xlabel('n')
hold on
plot(En)
legend('Diagonalization eigenvalues','Box exact energy values')
subplot(1,2,2)
plot(abs((En)'-e))
legend('divergence of eigenvalues')
axes('Position',[0.2 0.49 0.5 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0, 0.98,'Numeric and analytic eigenvalues ')
text(1, 0.98,'|E_{n}-e_{n}|')
xlabel('n')

%-------plotting numerical and analytical functions-----------
figure(2)
subplot(2,3,1)
plot(x,-u(:,1))
hold on
plot(x,sqrt(2/N)*cos(pi*x/dx/N))
xlabel('x')
legend('u_{1}','\psi_{1}')
subplot(2,3,2)
plot(x,u(:,2))
hold on
plot(x,sqrt(2/N)*sin(2*pi*x/dx/N))
xlabel('x')
legend('u_{2}','\psi_{2}')
subplot(2,3,3)
plot(x,-u(:,3))
hold on
plot(x,sqrt(2/N)*cos(3*pi*x/dx/N))
xlabel('x')
legend('u_{3}','\psi_{3}')
subplot(2,3,4)
plot(x,-u(:,4))
hold on
plot(x,sqrt(2/N)*sin(4*pi*x/dx/N))
xlabel('x')
legend('u_{4}','\psi_{4}')
subplot(2,3,5)
plot(x,-u(:,5))
hold on
plot(x,sqrt(2/N)*cos(5*pi*x/dx/N))
xlabel('x')
legend('u_{5}','\psi_{5}')
subplot(2,3,6)
plot(x,-u(:,6))
hold on
plot(x,sqrt(2/N)*sin(6*pi*x/dx/N))
xlabel('x')
legend('u_{6}','\psi_{6}')
axes('Position',[0.2 0.49 0.5 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.45, 0.98,'numerical u_{n} and analytical \psi_{n} solutions')

%-------plotting difference numerical and analytical eigenfunctions--------
figure(3)
subplot(2,3,1)
plot(x,abs(-u(:,1)'-sqrt(2/N)*cos(pi*x/dx/N)))
legend('|u_{1}-\psi_{1}|')
subplot(2,3,2)
plot(x,abs(-u(:,2)'+sqrt(2/N)*sin(2*pi*x/dx/N)))
legend('|u_{2}-\psi_{2}|')
subplot(2,3,3)
plot(x,abs(-u(:,3)'-sqrt(2/N)*cos(3*pi*x/dx/N)))
legend('|u_{3}-\psi_{3}|')
subplot(2,3,4)
plot(x,abs(-u(:,4)'-sqrt(2/N)*sin(4*pi*x/dx/N)))
legend('|u_{4}-\psi_{4}|')
subplot(2,3,5)
plot(x,abs(-u(:,5)'-sqrt(2/N)*cos(5*pi*x/dx/N)))
legend('|u_{5}-\psi_{5}|')
subplot(2,3,6)
plot(x,abs(-u(:,6)'-sqrt(2/N)*sin(6*pi*x/dx/N)))
legend('|u_{6}-\psi_{6}|')
axes('Position',[0.2 0.49 0.5 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.3, 0.98,'difference between numerical u_{n} and analytical \psi_{n} solutions |u_{n}-\psi_{n}|')

end

