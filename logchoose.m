function R = logchoose(n,i);
%This function computes the natural logarithm of "n choose i", or binomial coefficient 

A=0;
B=0;
C=0;

if i==0
    R=0;
else
    if (i<0) | (i>n)
        R = -5e200;
    else
        
        
        for j=1:n
            A=A+log(j);
        end
        
        for j=1:i
            B=B+log(j);
        end
        
        for j=1:n-i
            C=C+log(j);
        end
        
        R=A-B-C;
    end
end

