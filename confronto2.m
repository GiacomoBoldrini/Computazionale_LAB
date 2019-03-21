
k=(sin(1/2))^2;
K=ellipke(k);
T=4*K;
y0 = [1 0];

x1 = zeros(1,13);
err1 = zeros(1,13);
steps1 = zeros(1,13);
time1=zeros(1,13);
for j = 1:13
tol = 10^(-j);
x1(j) = 10^(-j);
opts = odeset('reltol',tol,'abstol',tol,'refine',1);
tic
[t,y] = ode45(@pendulum,[0 10*T],[1 0]',opts);
time1(j) = toc;
steps1(j) = length(t)-1;
err1(j) = max(abs(y(end,:)-y0));
end


subplot(3,2,1)
loglog(x,steps1,'--o');
hold on
subplot(3,2,2)
loglog(x,err1,'--o');
hold on
subplot(3,2,3)
loglog(x,time1,'--o');
hold on



x = zeros(1,13);
err = zeros(1,13);
steps = zeros(1,13);
time=zeros(1,13);

for j = 1:13
tol = 10^(-j);
x(j) = tol;
tic
F=@(q) -sin(q);
[t,q,p] = odesympvariable(F,[0 5*T],1,0,tol);
subplot(3,2,2)
loglog(t,0.5*p.^2-cos(q),'--+')
hold
time(j) = toc;
steps(j) = length(t)-1;
err(j) = max(abs(p(end)));
end

subplot(3,2,1)
loglog(x,steps,'--o');
legend('ode45','odesymp')
xlabel('tolerance')
ylabel('steps')
subplot(3,2,2)
loglog(x,err,'--o');
legend('ode45','odesymp')
xlabel('tolerance')
ylabel('error')
subplot(3,2,3)
loglog(x,time,'--o');
legend('ode45','odesymp')
xlabel('tolerance')
ylabel('time')
disp('As we can see from the plots odesym always requires more steps to compute the integration ')
disp('The error of odesymp is less than the one from ode45 up to 1.e-12 where we have a swap')
disp('The time required is higher for ode45 up to 1.e-8. Smaller tolerances requires odesymp more time than ode45')
function ydot = pendulum(t,y)
    ydot = [y(2); -sin(y(1))];
end
