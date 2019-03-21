function [dt,tout,qout,pout] = mysympint(Fa,tspan,q,p,dt,varargin)
% Symplectic midpoint
% the initial q and p must be column vectors.
t0 = tspan(1);
tfinal = tspan(2);
plotit = (nargout == 0);
t = t0;
x1 = 1/(2-2^(1/3));
x0 = 1 - 2*x1;
dtz=dt;

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
i=1;

% The main loop.
while i~=4
    if i==1 | i==3
        dt=dtz*x1;
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
            tout(end+1,1) = t; %#ok<*AGROW>
            qout(end+1,:) = q.';
            pout(end+1,:) = p.';
            end
        end

        if plotit
        plotter([],[],'done');
        end
        i=i+1;
    end
    if i==2
        dt=dtz*x0;
        while t < tfinal

             if 1.1*abs(dt) >= abs(tfinal - t)
                dt = tfinal - t;
             end
  
            t  = t + dt;
            q = q + 0.5*dt*p; %=q(n+1/2)
            a = Fa(q);%=a(q(n+1/2))
            p = p + dt*a;% p(n+1)=p(n)+dt*a(q(n+1/2))
            q = q + 0.5*dt*p; %q(n+1)=q(n)+dt*p(n)+0.5*dt^2*a(q(n+1/2))

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
        i=i+1; %non sto iterando l'operazione, lui mi stampa fuori solo l'ultima operazione fatat che equivale ad un passo di h*x1
    end
end
        