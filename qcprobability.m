function qcprobability(x,u,V,dx,e,N)
%qc probability shows the quantum-classical correspondence and it's called by 'finale'.
%As explained the semiclassical WKB approximation gives as a result a 
%probability to find the particle in a region between x and x+dx proportional
%to 1/p but it is wrong for low quantum numbers. Instead, for larger 
%quantum numbers, we get the classical interpretation satisfying the Bohr 
%correspondence principle.

k=1;
for j=[1, 5, 8 , 10]
i=find(V<=e(j));
subplot(2,2,k)
hold on
plot(x(i),(sqrt(2*e(j)/(2*(N*dx))))*dx*1./abs(sqrt(2*(e(j)-V(i)))),'--')
plot(x,abs(u(:,j)).^2)
legend('Classical','Quantum')
k=k+1;
axes('Position',[0.1 0.49 0.5 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Classical and quantum (|u_{j}|^{2}) probability distributions for j=1,5,8,10')
xlabel('x')
ylabel('P')
end

end