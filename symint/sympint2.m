function [tout,qout,pout] = sympint2(Fa,tspan,q,p,dt,varargin)
% Symplectic midpoint
% the initial q and p must be column vectors.

t0 = tspan(1);
tfinal = tspan(2);
plotit = (nargout == 0);
t = t0;

plotter = @odeplot;
%plotter = @odephas2;

% Initialize output.

if plotit
   plotter(tspan,[q;p;0.5],'init');
else
   tout = t;
   qout = q.';
   pout = p.';
end

% The main loop.

while t < tfinal

    if 1.1*abs(dt) >= abs(tfinal - t)
      dt = tfinal - t;
    end
  
    t  = t + dt;
    q = q + 0.5*dt*p;
    a = Fa(q);
	p = p + dt*a;
	q = q + 0.5*dt*p;

    if plotit
      plotter(t,[q;p;0.5],'');
    else
      tout(end+1,1) = t; %#ok<*AGROW>
      qout(end+1,:) = q.';
      pout(end+1,:) = p.';
   end
end

if plotit
   plotter([],[],'done');
end
