function [t,q,p] = sympint4(q,p,dt,fa)

	x1 = 1/(2-2^(1/3));
	x0 = 1 - 2*x1;

	
	[t,q,p] = sympint2(fa,q,p,x1*dt);
	[t,q,p] = sympint2(fa,q,p,x0*dt);
	[t,q,p] = sympint2(fa,q,p,x1*dt);
	

	
