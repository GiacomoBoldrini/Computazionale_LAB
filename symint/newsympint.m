function [tout3,qout3,pout3] = newsympint(Fa,tspan,q,p,dt,varargin)
% Symplectic midpoint
% the initial q and p must be column vectors.
x1 = 1/(2-2^(1/3));
x0 = 1 - 2*x1;
    
t0 = tspan(1);
tfinal = tspan(2);
plotit = (nargout == 0);
t = t0;

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
i=0;
while i~=4 
    case 1 
    dt=dt*x1;
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
      plotter(t,[q;p],'');
    else
      tout1(end+1,1) = t; %#ok<*AGROW>
      qout1(end+1,:) = q.';
      pout1(end+1,:) = p.';
    end
    end
    i=i+1;
    
    case 2
     dt=dt*x0;
    while tout1 < tfinal

    if 1.1*abs(dt) >= abs(tfinal - tout1)
      dt = tfinal - tout1;
    end
  
    tout1  = tout1 + dt;
    qout1 = qout1 + 0.5*dt*pout1;
    a = Fa(qout1);
	pout1 = pout1 + dt*a;
	qout1 = qout1 + 0.5*dt*pout1;

    if plotit
      plotter(tout1,[qout2;pout2],'');
    else
      tout2(end+1,1) = tout1; %#ok<*AGROW>
      qout2(end+1,:) = qout1.';
      pout2(end+1,:) = pout1.';
    end
    end
    i=i+1;
    
    case 3
     dt=dt*x1;
    while tout2 < tfinal

    if 1.1*abs(dt) >= abs(tfinal - tout2)
      dt = tfinal - tout2;
    end
  
    tout2  = tout2 + dt;
    qout2 = qout2 + 0.5*dt*pout2;
    a = Fa(q);
	pout2 = pout2 + dt*a;
	qout2 = qout2 + 0.5*dt*pout2;

    if plotit
      plotter(tout2,[qout2;pout2],'');
    else
      tout3(end+1,1) = tout2; %#ok<*AGROW>
      qout3(end+1,:) = qout2.';
      pout3(end+1,:) = pout2.';
    end
    end
    i=i+1;
    
    

end
end
