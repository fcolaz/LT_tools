%this script runs the finite length analysis of an LT code under
%inactivation decoding following [2]

% [2] Lázaro, Francisco, Gianluigi Liva, and Gerhard Bauch. "Inactivation 
% decoding of LT and Raptor codes: Analysis and code design." IEEE Trans.
% on Commun. 65.10 (2017): 4114-4127.

clear all
%close all





delta_v = (0:10:100);           % overhead values to be evaluated
x.k=100;                        % number of input symbols

%degree distribution initialization 
x.Omega= genSolitonDist('RSD', x.k, 0.02, 0.05); % ISD/RSD
%x.Omega = get_LT_dist_R10();   % degree distribution from R10 Raptor codes


x.max_degree = length(x.Omega); %maximum degree of the degree distribution


N_inact = zeros(size(delta_v)); % This vector will contain the number of 
                                % inactivations at each value of delta_v



% Now we precompute a look up table with the (logarithm of the) binomial
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
        
    N_inact(idx_delta) = get_n_inact(x);
    
end
toc

figure
plot(delta_v, N_inact)
hold on
%the next two lines are the result of simulating an RSD dist with k=100,
%c=0.02, and delta =0.05, collecting 10000 differente realizations
my_delta = 0:10:100;  
my_n_inact = [12.3739    6.8610    3.1306    1.3818    0.7063    0.3900    0.2382    0.1419    0.0990    0.0695    0.0452];
plot(my_delta,my_n_inact, '*')
grid minor
legend('analysis', 'simulation')
xlabel('\delta')
ylabel('number of inactivations')
title('number of inactivations vs \delta')




