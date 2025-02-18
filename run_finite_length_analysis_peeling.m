%this script implements the finite length analysis of an LT code under
%peeling decoding proposed in [1]

%[1] Karp, Richard, Michael Luby, and Amin Shokrollahi. "Finite length 
% analysis of LT codes." Proc of the Int. Symp. on Inf. Theory (ISIT) 2004. 


clear all
close all


delta_v = (0:10:100);           % overhead values to be evaluated
x.k=100;                        % number of input symbols

%degree distribution initialization 
x.Omega= genSolitonDist('RSD', x.k, 0.02, 0.05); % ISD/RSD
%x.Omega = get_LT_dist_R10();   % degree distribution from R10 Raptor codes


x.max_degree = length(x.Omega); %maximum degree of the degree distribution


P_f = zeros(size(delta_v));  % This vector will contain the probability 
%                              of decoding failure at each value of delta_v



% We precompute a look up table with the (logarithm of the) binomial
% coefficients.
m_max = round( x.k + delta_v(end));  
x.logchoose_matrix = get_logchoose_matrix(m_max); 


% We also precompute p_u, which depends only on the degree distribution
x. p_u = get_p_u(x);



% The following loop runs the finite length analysis for each value of
% delta

tic
for idx_delta = 1: length(delta_v)

    x.delta =delta_v(idx_delta);            
        
    [Pf_, Pf_u, ripple_u, cloud_u] = peeling_decoding_analysis(x);
    P_f(idx_delta)= Pf_;
    
end
toc

figure
semilogy(delta_v, P_f)
hold on
%the next two lines are the result of simulating an RSD dist with c=0.02,
%and delta =0.05, collecting 100 errors
my_delta = 0:10:100;  
my_p_f = [1, 0.996016, 0.838926, 0.543478, 0.316256, 0.215983, 0.14174, 0.0953107, 0.0711339, 0.0475694, 0.0371664];
plot(my_delta,my_p_f, '*')
grid minor
legend('analysis', 'simulations')
xlabel('\delta')
ylabel('P_f')
title('P_f vs \delta')

figure
plot([0:1:x.k], ripple_u, 'b')
hold on
plot([0:1:x.k], cloud_u, 'r')
xlabel('u')
ylabel('ripple/cloud size')
grid on
title('ripple/cloud size vs u')
legend('ripple(u)', 'cloud(u)', 'Location','northwest')


figure
plot([1:1:x.k], Pf_u(2:end), 'b')
xlabel('u')
ylabel('Prob. of decoder stopping at u')
grid on
title('P_f(u)')

