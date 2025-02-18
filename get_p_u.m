function p_u = get_p_u(x)
% Computes P_u according to proposition 3 in [2]
% [2] Lázaro, Francisco, Gianluigi Liva, and Gerhard Bauch. "Inactivation 
% decoding of LT and Raptor codes: Analysis and code design." IEEE Trans.
% on Commun. 65.10 (2017): 4114-4127.

% Params:  x:              Data structure describing an LT code  

%          x.k:            number of input symbols
%          x.Omega:        LT degree distribution
%          x.max_degree:   maximum degree (length) of x.Omega
%          x.logchoose_matrix:     look up table for the log of the 
%                                  binomial coefficient

% Output: p_u a vector of size (1 x (k+1))
%p_u(u+1) returns the probability that a randonly chosen output symbol
%leaves the cloud and enters the ripple in the transision from u to u-1
%active symbols (unprocessed symbols)

logchoose_matrix = x.logchoose_matrix;
k=x.k;
Omega = x.Omega;
p_u = zeros(1,k+1);


for idx= 3:k+1
    u=idx-1;
    numerator = 0;
    denominator_1 = 0;
    denominator_2 = 0;
    for d = 1:x.max_degree

        %we compute the numerator
        if (d>=2) && (d<=k-u+2)
            tmp = logchoose_matrix(k-u+1, d-2+1) - logchoose_matrix(k+1, d+1);
            numerator = numerator + Omega(d) * exp(tmp);
        end
        %we compute the 1st term of the denominator
        if d<= k-u+1
            tmp2 =  logchoose_matrix(k-u+1,d-1+1) - logchoose_matrix(k+1,d+1);
            denominator_1= denominator_1 + Omega(d) *  exp(tmp2);
        end
        
        %we compute the 2nd term of the denominator
        if d<= k-u
            tmp3 = logchoose_matrix(k-u+1,d+1) - logchoose_matrix(k+1,d+1);
            denominator_2 = denominator_2 + Omega(d) *  exp(tmp3);
        end
               
    end
    %we aggregate the results
    numerator = (u-1) * numerator;  
    denom = 1 - u *  denominator_1 - denominator_2;
    p_u(idx) = numerator/denom;    
end

% p_u is a probability, we ensure it takes a value between 0 and 1
p_u= min(p_u,1);
p_u=max(p_u,0);

