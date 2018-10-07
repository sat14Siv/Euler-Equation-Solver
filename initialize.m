function [rho, u, T, rhoWrite, uWrite, pWrite] = initialize(rhoL, rhoR, uL, uR, TL, TR, bc, iters, writeFreq, len)
    R = 287;
    rho = zeros(len-1,1);
    rho(1: 0.5*(size(rho,1))) = rhoL;        
    rho(0.5*(size(rho,1)+2) : end) = rhoR;

    u = zeros(len-1,1);
    u(1: 0.5*(size(rho,1))) = uL;        
    u(0.5*(size(rho,1)+2) : end) = uR;
    u(1) = bc*uL;
    u(end) = bc*uR;
    
    T = zeros(len-1,1);
    T(1: 0.5*(size(T,1))) = TL;        
    T(0.5*(size(T,1)+2) : end) = TR;
    
    % Write variables
    rhoWrite = zeros(len-2, iters/writeFreq + 1);
    rhoWrite(1,1) = 0 ; rhoWrite(2:end, 1) = rho(2:end-1);
    uWrite = zeros(len-2, iters/writeFreq + 1);
    uWrite(1,1) = 0 ; uWrite(2:end, 1) = u(2:end-1);
    pWrite = zeros(len-2, iters/writeFreq + 1);
    pWrite(1,1) = 0 ; pWrite(2:end, 1) = rho(2:end-1).*T(2:end-1)*R;