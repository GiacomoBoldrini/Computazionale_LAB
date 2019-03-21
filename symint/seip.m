function [q,p]=seip(dqdt,dpdt,tspan,q0,p0,varargin)
%SEIP Symplectic Euler solver, implicit method.
%  [Q P] = SEIQ(DQDT,DPDT,TSPAN,Q0,P0) solves Hamilton's equations dq/dt = 
%  dH/dp, dp/t = - dH/dq with H = H(q,p) using an implicit Euler solver by 
%  stepping from T0 to T1 to T2, etc. Functions DQDT(Q,P) and DPDT(Q,P) must
%  return dH/dp and - dH/dq in the form of a N-dimensional column vectors.
%  Vectors Q0 and P0 are the initial conditions at T0. Each row in the 
%  solution arrays Q and P corresponds to a time specified in TSPAN.  
%
%  [Q P] = SEIP(DQDT,DPDT,TSPAN,Q0,P0,VARARGIN) passes the additional
%  parameters R1, R2, ... to functions DQDT(Q,P,R1,R2,...) and DPDT(Q,P,R1,
%  R2,...).
%
%  This numerical integrator treats the p-variable by the implicit Euler 
%  method and the q-variable by the explicit Euler method.
%  
%  For separable Hamiltonians, i.e H(q,p) = T(p) + V(q), the method is
%  explicit. Use SEEP for efficiency.
%
%  See also SEEQ, SEEP, SEIQ, GLS.

%  Francisco J. Beron-Vera, 28-Mar-2005
%  $Revision: 1.0 $  $Date: 28-Mar-2005 14:58:31 $

N = length(q0);  
Nt = length(tspan); 
hs = diff(tspan);
q = zeros(N,Nt); q(:,1)=q0;
p = zeros(N,Nt); p(:,1)=p0;
for nt = 2:Nt
   t = tspan(nt-1); if ~mod(nt,100), disp(['t = ' num2str(t)]), end
   h = hs(nt-1);
	for n = 1:N
		p(n,nt) = fzero(@(g) ...
         g - p(n,nt-1) - h * feval(dpdt, q(n,nt-1), g, varargin{:}), ...
         p(n,nt-1));
   end
	q(:,nt) = q(:,nt-1) + h * feval(dqdt, q(:,nt-1), p(:,nt), varargin{:}); 
end
q = q.';
p = p.';