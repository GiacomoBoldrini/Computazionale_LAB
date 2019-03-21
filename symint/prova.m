function    [probab,alias]=prova(p)
l=[1 2];
s=[3 4];
probab= [];
alias = [];
p=p*4

% while isempty(l)== 0 || isempty(s) == 0
%     l(1)=[];
%     s(1)=[];
%     i=i+1;
% end
 while isempty(s) == 0 & isempty(l) == 0
        z = s(1)
        k = l(1)
        s(1) = []
        l(1) = []
        probab(z) = p(z)
        alias(z) = k
        p(k) = p(k) + p(z)-1 
        if p(k) < 1
            s = [k,s]
        else
            l = [k,l]
        end
 end
end