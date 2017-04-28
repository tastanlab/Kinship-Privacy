function [x,fval, exitflag] = ga_lp2()

phi = dlmread('phi.txt','',0,0) #the maximum phi value found among family members before adding the newest member
M = dlmread('kin_constraints.txt','', 0,0 )  ;
s = (size(M,2)-2)/2;

x0 = [ zeros(s+1,1)];
lb = [zeros(s,1); phi];
size_constraints = dlmread('size_constraints.txt');
ub =[size_constraints ; 0.5];
IntCon=[1:s];

Aeq = [];
beq = [];
A = [];
b = [];

if(size(dir('*.vcf'),1) > 2)
    [x,fval,exitflag,output]=ga(@myobj,s+1,A,b,Aeq,beq,lb,ub,@myconstr, IntCon);   
else
    [x,fval,exitflag,output]=ga(@myobj2,s+1,A,b,Aeq,beq,lb,ub,@myconstr2, IntCon); 
end
fval
fileID = fopen('phi.txt','wt');
fprintf(fileID, '%f', fval');
fclose(fileID);

function f = myobj(x)
M = dlmread('kin_constraints.txt','', 1,0 )  ;
s = (size(M,2)-2)/2;
f=x(s+1);
end
%---------------------------------------------------------
function [c, ceq] = myconstr(x)
M = dlmread('kin_constraints.txt','', 0,0 ) ;
    s = (size(M,2)-2)/2;
    r = size(M,1);
    c = [];
    %% every constraint
    for i = 1:r
        a = 0;
        for j = 1:s
            a = a + x(j) * M(i,j);
            a = a + x(j) * M(i,j+s) * x(s+1);
        end
        a = a + x(s+1) * M(i,2*s+2);
        a = a - M(i,(2 * s + 1));
        c = [c; -a];
    end
    
    ceq=[];
end  
%---------------------------------------------------------
function f = myobj2(x)
f=x(2);
end
%---------------------------------------------------------    
function [c, ceq] = myconstr2(x)
c = [ -M(1)*x(1) - M(2)*x(1)*x(2) + M(4)*x(2) + M(3)];
ceq = [];
end
end
