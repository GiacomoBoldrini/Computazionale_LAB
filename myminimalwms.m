%% init simulation space
N = 2048;
%N = 1023;  % 2^m-1 to use dst() below
mybox = [-45 45];
x = linspace(mybox(1),mybox(2),N)';
dx = x(2)-x(1); 
g = squarelattice(N);
L = g.laplacian/dx/dx;
% wavenumbers and khat's
n = floor(N/2);
nn = floor((N-1)/2);
k = (2*pi/(dx*N))*(-n:nn)';
khat = (2/dx)*sin((pi/N)*(0:N-1)');
%% init vars and plot
%V = 0.25*x.^4;
%V = 0.5*x.^2 + 0.005*x.^4;
%V = zeros(N,1);
%V = 15*exp(-(x).^2/100);
%V = 10*exp(-(x).^2/100);
V=15*(-heaviside(x)+heaviside(x+3));
k0 = 5;
%psi = exp(-(x-15).^2/2);
psi = exp(1i*k0*x).*exp(-(x+35).^2/8);
psi = psi/(sqrt(dx)*norm(psi));
%kpsi = sfft(psi);
%clf
% s = 0.005;
s = 0.05;
figure(1)
clf
plot(x,s*V)
axis([mybox(1) mybox(2) -1 1.25])
hold on
box on
h = plot(x,real(psi));
% %figure(2)
% clf
% %hk = plot(k,real(kpsi));
% axis([mybox(1) mybox(2) -50 80])
% hold on
% box on
%kpsi = sfft(psi);
%h = plot(x,s*real(kpsi));
%axis([-1 1 -1 1])
grid
nplot = 5;

%% init run
 kk = fftshift(k);
% kk = khat;
% kk=kd; % non devo fare lo fftshift perchè il sinft prende come indici [-1;N]
% esattamente come è stato creato il nostro indice
%kk = khatd;
dt = 0.01;
%dt = -1i*0.01;
nt = 2300;
upsi = psi;
%upsi = conj(upsi);

%% run
for j = 1:nt
    %upsi = upsi - (dx*psi0'*upsi)*psi0;
    upsi = exp(-1i*dt*V/2).*upsi;
    upsi = ifft(exp(-1i*dt*kk.^2/2).*fft(upsi));
    %upsi = idst(exp(-1i*dt*kk.^2/2).*dst(upsi));
    %upsi = sinft(exp(-1i*dt*kk.^2/2).*sinft(upsi)); %la trasformata di 
    % Fourier con i seni è ortogonale: coincide quindi con l'inversa 
    % e non ho quindi la "invsinft"
    upsi = exp(-1i*dt*V/2).*upsi;
    %plot(k,kpsi)
    upsi = upsi/(sqrt(dx)*norm(upsi));
    if mod(j,nplot) == 0
        set(h,'Ydata',abs(upsi));
        %set(h,'Ydata',real(upsi));
        %kpsi = sfft(upsi);
        %set(hk,'Ydata',abs(kpsi).^2/1000);
        drawnow
        %norm(upsi)
        %pause(0.2)
    end
end
%% average kinetic energy

kpsi = sfft(upsi);
dk = 2*pi/dx;
T = sum(k.^2.*abs(kpsi).^2)/sum(abs(kpsi).^2)
%% trasmission coefficient
Tr = sum(abs(upsi(x>0)).^2)/sum(abs(upsi).^2) 
%% reflection coefficient
R = sum(abs(upsi(x<0)).^2)/sum(abs(upsi).^2) 