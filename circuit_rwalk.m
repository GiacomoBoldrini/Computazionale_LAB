function s = circuit_rwalk(rvals)

if nargin == 0
    rvals = rand(1,8);
end

%% prep
B = [ 1 -1  0  0  0 
      1  0 -1  0  0 
     -1  0  0  1  0 
      0 -1  1  0  0 
      0  0 -1  1  0 
      0  1  0  0 -1 
      0  0  1  0 -1 
      0  0  0 -1  1];
  
syms vs is r12 r13 r14 r23 r34 r25 r35 r45
r = [r12 r13 r14 r23 r34 r25 r35 r45];
g = 1./r;
G = B'*diag(g)*B; % conductivity matrix
I = [0 0 1 0 -1]';
v = G\I;  % voltage (same as v = [G(1:4,1:4)\i0(1:4) ; 0] in this special case)

%% Let us check
disp('-----------------------------')
disp('')
disp(['v(5) = ' char(v(5))])
I0 = simplify(G*[1 1 1 1 1]')
I = simplify(G*v)
% assign to output struct
s.I0 = I0;
s.I = I;
disp('-----------------------------')
disp('')
disp('---- press any key to continue ---- ')
pause

%% stochastic matrices and (di)graphs
A = diag(diag(G)) - G;
g0 = graph(double(subs(A,r,ones(1,8)))');
x = [.5 .2 .5 .8 .7];
y = [.6 .6 .9 .7 .2];
figure(1)
lbl = {'r12' 'r13' 'r14' 'r23' 'r34' 'r25' 'r35' 'r45'};
plot(g0,'XData',x,'YData',y,'EdgeLabel',lbl)
disp('figure 1: graph of conductances')
disp('---- press any key to continue ---- ')
pause

%%
%W0 = A./(ones(5,1)*diag(G).');
W0 = A./diag(G).';
%disp(['sum(W0) = ' char(simplify(sum(W0)))])
dg0 = digraph(double(subs(W0,r,ones(1,8)))');
figure(2)
lbl = {'w21' 'w31' 'w41' 'w12' 'w32' 'w52' 'w13' 'w23' ...
       'w43' 'w53' 'w14' 'w34' 'w54' 'w25' 'w35' 'w45'};
plot(dg0,'XData',x,'YData',y,'EdgeLabel',lbl,'ArrowSize',8)
disp('figure 2: digraph without generator')
disp('---- press any key to continue ---- ')
pause

%%
W = W0; W(:,5) = [0 0 1 0 0]';
dg = digraph(double(subs(W,r,ones(1,8)))');
figure(3)
lbl = {'w21' 'w31' 'w41' 'w12' 'w32' 'w52' 'w13' 'w23' ...
       'w43' 'w53' 'w14' 'w34' 'w54' '1'};
plot(dg,'XData',x,'YData',y,'EdgeLabel',lbl,'ArrowSize',8)
disp('figure 3: digraph with generator')
disp('---- press any key to continue ---- ')
pause

%% check with input resistences
W0 = double(subs(W0,r,rvals));
[ev,e] = eig(W0);
[~,I] = sort(diag(e),'descend'); ev = ev(:,I);
u0 = ev(:,1)/sum(ev(:,1));
double([u0 subs(diag(G)/sum(diag(G)),r,rvals)])    
%%
disp('-----------------------------')
disp('')
W = double(subs(W,r,rvals));
[ev,e] = eig(W);
[~,I] = sort(diag(e),'descend'); ev = ev(:,I);
u = ev(:,1)/sum(ev(:,1))

%%
disp('-----------------------------')
disp('')
d = double(subs(diag(G),r,rvals));
% potential from Ohm.s law v = G\I; 
vvals(:,1) = double(subs(v,r,rvals));
% potential from random walk equilibrium
vvals(:,2) = [u(1:4,1)./(u(5,1)*d(1:4)); 0];
% comparison
vvals

% assign to output struct
s.W0 = W0;
s.W = W;
s.v = vvals;



