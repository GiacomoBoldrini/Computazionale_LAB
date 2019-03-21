function  [tout1,qout1,pout1,tout,qout,pout]=mysympint2(Fa,tspan,q,p,dt,varargin)
% Symplectic midpoint
% the initial q and p must be column vectors.
q01 = q;
p01 = p;
t0 = tspan(1);
tfinal = tspan(2);
plotit = (nargout == 0);
t = t0;
c1=1/(2*(2-2.^(1/3)));
c2=(1-2.^(1/3))/(2*(2-2.^(1/3)));
c3=c2;
c4=c1;

d1=1/(2-2.^(1/3));
d2=-(2.^(1/3))/(2-2.^(1/3));
d3=d1;
d4=0;

rtol = 1.e-3;
atol = 1.e-6;
tdir = sign(tfinal - t0);
threshold = atol / rtol;
s1 = Fa(q, varargin{:});
r = norm(s1./max(abs(p),threshold),inf) + realmin;
h = tdir*0.8*rtol^(1/4)/r;
hmax = abs(0.1*(tfinal-t0));
q02 = q;
p02 = p;
q=0;
p=0;
q1=0;
p1=0;


%plotter = @odeplot;
plotter = @odephas2;

% Initialize output.

if plotit
   plotter(tspan,[q;p],'init');
else
   tout = t;
   qout = q.';
   pout = p.';
   
   tout1 = t;
   qout1 = q1.';
   pout1 = p1.';
end
% The main loop.
while t < tfinal

   hmin = 16*eps*abs(t);
   if abs(dt) > hmax, dt = tdir*hmax; end
   if abs(dt) < hmin, dt = tdir*hmin; end
   
   if 1.1*abs(dt) >= abs(tfinal - t)
   dt = tfinal - t;
   end
  
    %attempt a step
    
    tnew  = t + dt;
    q = q01 + c1*dt*p01;
    a = Fa(q);
	p = p + d1*dt*a;
	q = q + c2*dt*p;
    a = Fa(q);
	p = p + d2*dt*a;
	q = q + c3*dt*p;
    a = Fa(q);
	p = p + d3*dt*a;
	q = q + c4*dt*p;
    
    a = Fa(q02);
    p1 = p02 + c1*dt*a;
	q1 = q02 + d1*dt*p;
    a = Fa(q1);
	p1 = p1 + c2*dt*a;
	q1 = q1 + d2*dt*p1;
    a = Fa(q1);
	p1 = p1 + c3*dt*a;
	q1 = q1 + d3*dt*p1;
    a = Fa(q1);
    p1 = p1 + c4*dt*a;
    q1 = q1 + d4*dt*p1;
    
    err = max(norm(q-q1,inf),norm(p-p1,inf)) + realmin;
% Accept the solution if the estimated error is less than the tolerance
     if err <= rtol
      t = tnew;
      q01 = q;
      q02 = q1;
      if plotit
         if plotfun(q,p,'');
            break
         end
      else
         tout(end+1,1) = tnew;
         qout(end+1,:) = q.';
         pout(end+1,:) = p.';
         
         tout1(end+1,1) = t;
         qout1(end+1,:) = q1.';
         pout1(end+1,:) = p1.';
      end
     end
     dt = dt*min(5,0.8*(rtol/err)^(1/4)); %passo variabile
     if abs(h) <= hmin %se il passo è minore dello step minimo
      warning('Step size %e too small at t = %e.\n',dt,t);
      t = tfinal;
     end    
         
     
end 
%     if plotit
%       plotter(t,[q;p],'');
%     else
%       tout(end+1,1) = t; %#ok<*AGROW>
%       qout(end+1,:) = q.';
%       pout(end+1,:) = p.';
%       
%       tout1(end+1,1) = t;
%       qout1(end+1,:) = q1.';
%       pout1(end+1,:) = p1.';
%    end
% end

if plotit
   plotter([],[],'done');
   
end
