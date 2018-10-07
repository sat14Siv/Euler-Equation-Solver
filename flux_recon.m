function [R1, R2, R3] = flux_recon(M, ho, rho, p, u, imn, imx, a, fl_rec)
    % Flux Formulation Parameters
    alpha_p = 0.5*(1 + sign(M));
    alpha_n = 0.5*(1 - sign(M));
    beta = -max(0, 1 - floor(abs(M)));
    D_p = 0.25*(1 + M).^(2).*(2 - M);
    D_n = 0.25*(M - 1).^(2).*(2 + M);
    Di = alpha_p.*(1 + beta) - beta.*D_p;
    Di1 = alpha_n.*(1 + beta) - beta.*D_n;
    Cvl_p = alpha_p.*(1 + beta).*M - 0.25*beta.*((1 + M).^2);
    Cvl_n = alpha_n.*(1 + beta).*M + 0.25*beta.*((1 - M).^2);            
    
    Cvl_p1 = max(0,(Cvl_p(1:imx) + Cvl_n(2:imx+1)));   % AUSM
    Cvl_n1(2:imx+1) = min(0,(Cvl_p(1:imx) + Cvl_n(2:imx+1)));
    Cvl_p1(imx+1)=0;Cvl_n1(1) = 0;
    
    if fl_rec == 1
        Cvl_p = Cvl_p1;
        Cvl_n = Cvl_n1;
    end
    
    fl_q1 = zeros(size(M,1)-1,1);fl_q2 = zeros(size(M,1)-1,1);
    fl_q3 = zeros(size(M,1)-1,1);   
    R1 = zeros(size(M,1)-2,1);R2 = zeros(size(M,1)-2,1);
    R3 = zeros(size(M,1)-2,1);
    % Flux
    for i = imn-1: imx
        a_av = 0.5*(a(i) + a(i+1));
        fl_q1(i) = rho(i)*a_av*Cvl_p(i)*1  + rho(i+1)*a_av*Cvl_n(i+1)*1;
        fl_q2(i) = rho(i)*a_av*Cvl_p(i)*u(i)  + rho(i+1)*a_av*Cvl_n(i+1)*u(i+1)...
            + (Di(i)*p(i)+Di1(i+1)*p(i+1))*1;
        fl_q3(i) = rho(i)*a_av*Cvl_p(i)*ho(i)  + rho(i+1)*a_av*Cvl_n(i+1)*ho(i+1);
    end
    % Residue
    for i = imn: imx
        R1(i-1) = fl_q1(i) - fl_q1(i-1);
        R2(i-1) = fl_q2(i) - fl_q2(i-1);
        R3(i-1) = fl_q3(i) - fl_q3(i-1);
    end