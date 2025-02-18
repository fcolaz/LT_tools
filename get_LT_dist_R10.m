function Omega = get_LT_dist_R10()
% This function returns the LT degree distribution used in R10 Raptor
% codes as defined in [3]
% [3 ]Luby, M., et al. "RFC 5053: Raptor forward error correction scheme 
% for object delivery." (2007)

Omega= zeros(1,40);
Omega(1) = 0.009766579;
Omega(2) = 0.459042549;
Omega(3) = 0.210964203;
Omega(4) = 0.11339283;
Omega(10) = 0.11134243;
Omega(11) = 0.079863548;

Omega(end) = 1- sum(Omega(1:end-1));








