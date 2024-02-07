function [rho, omegam, omegap, mudynamic, sigk, omegak,...
    omegaeps, cep1, cep2, cmu, sigeps, kappa, B, E] = constantes(Re)

    rho = 1;    
    omegam = 0.5;   
    omegap = 0.3;   
    mudynamic = 1/Re;
    sigk = 1;   
    omegak = 0.3;   
    omegaeps = 0.3;
    cep1 = 1.44;    
    cep2 = 1.92;    
    cmu = 0.09;
    sigeps = 1.3;   
    kappa = 0.41;   
    B = 5.5;

    E = exp(kappa*B);           % log-law constant 
    
end