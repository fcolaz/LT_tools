function [ distribution ] = genSolitonDist( type, k, c, delta )
% Generates either a ideal soliton distribution or a robust soliton
% distribution as defined in 
% Luby, Michael. "LT codes." Proc of the IEEE Symp. on Found. of Comp, Science, 2002.
%
% Parameters:   type:       Defines the type of the soliton distribution.
%                           'ISD' for ideal soliton distribution or 'RSD'
%                           for robust soliton distribution is possible.
%               k:          number of input symbols
%               c:          RSD parameter c
%               delta:      RSD parameter delta

if strcmp(type, 'ISD')
    distribution = zeros(1,k);
    distribution(1) = 1/k;
    for d = 2:k                       % generate ideal soliton distribution
        distribution(d) = 1/(d*(d-1));% according to appropriate equation
    end    
end

if strcmp(type, 'RSD')
    ISD = zeros(1,k);                 
    ISD(1) = 1/k;
    for d = 2:k
        ISD(d) = 1/(d*(d-1));
    end

    tau = zeros(1,k);
    S = c*log(k/delta)*sqrt(k);  

    for d = 1:round(k/S)-1       
        tau(d) = S/(k*d);        
    end
    tau(round(k/S)) = (S/k)*log(S/delta);
    for d = round(k/S)+1:k
        tau(d) = 0;
    end
    Z = sum(ISD + tau);
    distribution = (ISD + tau)/Z;
    
end

if strcmp(type, 'RSDtrunc')
    ISD = zeros(1,k);            
    ISD(1) = 1/k;
    for d = 2:k
        ISD(d) = 1/(d*(d-1));
    end

    tau = zeros(1,k);
    S = c*log(k/delta)*sqrt(k);  

    for d = 1:round(k/S)-1       
        tau(d) = S/(k*d);        
    end
    tau(round(k/S)) = (S/k)*log(S/delta);
    for d = round(k/S)+1:k
        tau(d) = 0;
    end
    Z = sum(ISD + tau);
    distribution = (ISD + tau)/Z;

    distribution = distribution(1:round(k/S));
    distribution = distribution./(sum(distribution));
  
end

end

