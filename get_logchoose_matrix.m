function logchoose_matrix = get_logchoose_matrix(m)
% This function returns a look-up table of the binomial coefficient.
% In particular the look-up talbe has dimensions (m+1) x (m+1) and 
% logchoose_matrix (n,i) corresponds to "(n-1) choose (i-1)"



logchoose_matrix = -1e6 * ones(m+1,m+1);
for i=0:m    
        logchoose_matrix (i+1,1) = 0;    
    for j=1:i
        logchoose_matrix(i+1,j+1) = logchoose_matrix(i+1,j)  + log(i-j+1) - log(j);        
    end
end
        
   
end