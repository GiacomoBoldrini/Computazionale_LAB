function [tout,qout,pout] = mysympint3(Fa,tspan,q,p,varargin)
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
d4=0;

rtol = 1.e-3;
atol = 1.e-6;
tdir = sign(tfinal - t0);
threshold = atol / rtol;
dtmax = abs(0.1*(tfinal-t0));

%initial step size
s1 = Fa(q);
r = norm(s1./max(max(abs(q),abs(p)),threshold),inf) + realmin;
dt =0.1; %tdir*0.8*rtol^(1/4)/r;

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
q0=q;
p0=p;

% The main loop.
while t < tfinal
    %step size
    dtmin = 16*eps*abs(t);
   if abs(dt) > dtmax, dt = tdir*dtmax; end
   if abs(dt) < dtmin, dt = tdir*dtmin; end
   
   %STEP FINALE
    if 1.1*abs(dt) >= abs(tfinal - t)
      dt = tfinal - t;
    end
    q0=q;
    p0=p;
    i=0;
    %SIMPLECTIC EULER 1
      
    t  = t + dt;

    if i==0 
    q=q0;
    p=p0;
    end

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
    
    q_SE1=q;
    p_SE1=p;
    
    
    %SIMPLECTIC EULER 2
    
    if i==0 
    q=q0;
    p=p0;
    end
    
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
    
    q_SE2=q;
    p_SE2=p;

 %BISOGNA CALCOLARE L'ERRORE
  
   err = max(norm(q_SE2-q_SE1,inf),norm(p_SE2-p_SE1,inf)) + realmin;
 % Accept the solution if the estimated error is less than the tolerance.

   if err <= rtol
      qnew = (q_SE1+q_SE2)/2;
      pnew = (p_SE1+p_SE2)/2;
      tnew=t;
      if plotit
         if plotfun(tnew,qnew,'');
            break
         end
      else
         tout(end+1,1) = tnew;
         qout(end+1,:) = qnew.';
         pout(end+1,:) = pnew.';
      end
      
   end
   
   % Compute a new step size.

   dt =0.1; %dt*min(5,0.8*(rtol/err)^(1/4));    
    
   
   


         
         i=i+1;
         
  
end

if plotit
   plotter([],[],'done');
end