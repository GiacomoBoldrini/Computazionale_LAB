function  [tout,qout,pout]=simpintt(Fa,tspan,q,p,dt,varargin)
% Symplectic midpoint
% the initial q and p must be column vectors.
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

rtol = 1.e-3;
atol = 1.e-6;
tdir = sign(tfinal - t0);
threshold = atol / rtol;
s1 = Fa(q, varargin{:});
r = norm(s1./max(abs(p),threshold),inf) + realmin;
h = tdir*0.8*rtol^(1/4)/r;

%plotter = @odeplot;
plotter = @odephas2;

% Initialize output.

if plotit
   plotter(tspan,[q;p],'init');
else
   tout = t;
   qout = q.';
   pout = p.';
end
% The main loop.
while t < tfinal

   hmin = 16*eps*abs(t);
   if abs(dt) > hmax, dt = tdir*hmax; end
   if abs(dt) < hmin, dt = tdir*hmin; end
   
   if 1.1*abs(dt) >= abs(tfinal - t)
   dt = tfinal - t;
    end
  
    t  = t + dt;
    q = q + c1*dt*p;
    a = Fa(q);
	p = p + d1*dt*a;
	q = q + c2*dt*p;
    a = Fa(q);
	p = p + d2*dt*a;
	q = q + c3*dt*p;
    a = Fa(q);
	p = p + d3*dt*a;
	q = q + c4*dt*p;
    
    pparg = p;
    qparg = q;
    
    
    
    a = Fa(q);
    p = p + c1*dt*a;
	q = q + d1*dt*p;
    a = Fa(q);
	p = p + c2*dt*a;
	q = q + d2*dt*p;
    a = Fa(q);
	p = p + c3*dt*a;
	q = q + d3*dt*p;
    a = Fa(q);
    p = p + c4*dt*a;
    q = q + d4*dt*p;

    if plotit
      plotter(t,[q;p],'');
    else
      tout(end+1,1) = t; %#ok<*AGROW>
      qout(end+1,:) = q.';
      pout(end+1,:) = p.';
   end
end

if plotit
   plotter([],[],'done');
   
end
