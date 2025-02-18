function [ n_inact ] = get_n_inact( x )
% This function takes an input a data structure "x" defining the parameters
% of an LT code and runs the finite length analysis of the code under
% inactivation decoding. The output is the expected number of inactivations
% "n_inact"
% Params:  x:              Data structure describing an LT code  

%          x.k:            number of input symbols
%          x.Omega:        LT degree distribution
%          x.max_degree:   maximum degree (length) of x.Omega
%          x.logchoose_matrix:     look up table for the log of the 
%                                  binomial coefficient
%          x.p_u:          p_u according to proposition 3 in [2]

% [2]   Lázaro, Francisco, Gianluigi Liva, and Gerhard Bauch. "Inactivation 
%       decoding of LT and Raptor codes: Analysis and code design." IEEE Trans.
%       on Commun. 65.10 (2017): 4114-4127.


%parameters to speed up the analysis

k_std = 5;% used to discard values of a and b which are very unlikely.
% only values of a and b which are less than k_std standard
% deviations away from the mean will be evaluated

thres = 1e-12; % decoder states whose probability is lower than thres
% will not be evaluated

r_max = round(x.k+x.delta);  %Maximum ripple size considered



k=x.k;
p_u=x.p_u;
logchoose_matrix = x.logchoose_matrix;


delta = x.delta;




m=k+delta;

P_c_r_u_prev = zeros( m+1, r_max +1 );

%here we compute the initial state of the decoder
o1= x.Omega(1);
for r = 0:r_max
    c = m-r;    
    tmp = logchoose_matrix(m+1,c+1) + c * log(1-o1) + (m-c) * log(o1);
    P_c_r_u_prev( c+1 , r +1) = exp(tmp);
end




N_inact_j = zeros (1,k+1);


P_c_r_u_next = P_c_r_u_prev;
N_inact_j(k+1) = sum(P_c_r_u_next(:,1));

for u= k:-1:1 % u is the number of unresolved symbols
    P_c_r_u_prev = P_c_r_u_next;
    P_c_r_u_next = zeros( m+1, r_max +1 );
    for r = 0: r_max %previous ripple
        %a --> a+1 symbols leave the ripple and become of reduced degree 0
        % a is distributed as Binomial (r, 1/u)
        % Next, we compute the minimum and maximum value of a that we will
        % consider
        a_mean = r * 1/u;
        a_var = r * (1/u) * (1-1/u);
        a_std = sqrt(a_var);

        min_a = a_mean - k_std * a_std;
        max_a = a_mean + k_std * a_std;
        min_a = max (0, floor(min_a));
        max_a = min (r-1, ceil(max_a));
        max_a = max(max_a, 0);
        min_a =0;
        for c = 0:m-r %previous cloud
            %b --> number of symbols leave the cloud enter the ripple. 
            % b is distributed as Binomial (c,p_u(u+1));
            % Next, we compute the minimum and maximum value of b that we 
            % will  consider
            b_mean = c * p_u(u+1);
            b_var = c * p_u(u+1) * (1-p_u(u+1));
            b_var = max(b_var, 0);
            b_std = sqrt(b_var);
            min_b = b_mean - k_std * b_std;
            max_b = b_mean + k_std * b_std;
            min_b = max (0, floor(min_b));
            max_b = min (c, ceil(max_b));
            if (P_c_r_u_prev(c+1,r+1) > thres) %if the state probability is too low, we skip it
                if r>0 
                    for a=min_a:max_a
                        tmp_a = logchoose_matrix(r-1+1,a+1) + a* log(1/(u)) + (r-1-a)* log(1-1/(u)+1e-200);
                        b1= -(r-a-1);
                        min_b_=max(min_b, b1);
                        b2 = r_max -r +a +1;
                        max_b_=min(max_b, b2);
                        for b = min_b_:max_b_
                            % a --> a+1 symbols leave the ripple
                            % b --> b symbols leave the cloud and enter the ripple

                            r_next = r -a-1 +b; % next ripple
                            c_next = c -b;  % next cloud
                            tmp = logchoose_matrix(c+1,b+1) +b * log( p_u(u+1) +1e-200) +(c-b) * log(1-p_u(u+1)+1e-200) + tmp_a;
                            P_c_r_u_next(c_next+1,r_next+1) = P_c_r_u_next(c_next+1,r_next+1) + exp( tmp) * P_c_r_u_prev(c+1,r+1);

                        end %for b
                    end %for a
                else
                    % since r=0, we have a= 0, and we only loop over b
                    for b = min_b:max_b
                        r_next = r +b;
                        c_next = c -b;
                        if ((r_next < r_max) && (c_next >=0))
                            tmp = logchoose_matrix(c+1,b+1) +b * log( p_u(u+1) +1e-30) +(c-b) * log(1-p_u(u+1)+1e-30) ;
                            P_c_r_u_next(c_next+1,r_next+1) = P_c_r_u_next(c_next+1,r_next+1) + exp( tmp) * P_c_r_u_prev(c+1,r+1);
                        end
                    end
                end
            end
        end
    end

    N_inact_j(u) = sum(P_c_r_u_next(:,1));   

end

n_inact = sum(N_inact_j(2:end));

end