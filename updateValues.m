function [rho, ho, T, a, M, p, u] = updateValues(q1_upd, q2_upd, q3_upd)
    R = 287; gamma = 1.4;
    rho = q1_upd; u = q2_upd./q1_upd; p = (q3_upd - 0.5*(q2_upd.^2)./q1_upd)*(gamma-1);  
    ho = q3_upd./q1_upd + p./rho;
    T = p./(rho*R);
    a = (gamma*R*T).^(0.5);
    M = u./a; 
    