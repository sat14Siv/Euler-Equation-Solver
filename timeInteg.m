function [q1_upd, q2_upd, q3_upd] = timeInteg(tg, del_x, q1, q2, q3, R1,...
    R2, R3, imn, imx, time_integ, fl_rec, bc)
    q1_upd = zeros(size(q1));q2_upd = zeros(size(q1));
    q3_upd = zeros(size(q1));
    
    if (time_integ == 0)
        q1_upd(imn:imx) = q1(imn:imx) - (tg/del_x)*R1;
        q2_upd(imn:imx) = q2(imn:imx) - (tg/del_x)*R2;
        q3_upd(imn:imx) = q3(imn:imx) - (tg/del_x)*R3;       
    end
    
    if (time_integ == 1)    % RK4
       R11 = R1; R12 = R2; R13 = R3;
       k11 = -R11/del_x; k12 = -R12/del_x; k13 = -R13/del_x;
       q1_upd(imn:imx) = q1(imn:imx) + 0.5*tg*k11;
       q2_upd(imn:imx) = q2(imn:imx) + 0.5*tg*k12;
       q3_upd(imn:imx) = q3(imn:imx) + 0.5*tg*k13;
       [q1_upd, q2_upd, q3_upd] = updateBC(q1_upd, q2_upd, q3_upd, imn, imx, bc);
       [rho, ho, ~, a, M, p, u] = updateValues(q1_upd, q2_upd, q3_upd);
       
       [R21, R22, R23] = flux_recon(M, ho, rho, p, u, imn, imx, a, fl_rec);
       k21 = -R21/del_x; k22 = -R22/del_x; k23 = -R23/del_x;
       q1_upd(imn:imx) = q1(imn:imx) + 0.5*tg*k21;
       q2_upd(imn:imx) = q2(imn:imx) + 0.5*tg*k22;
       q3_upd(imn:imx) = q3(imn:imx) + 0.5*tg*k23;
       [q1_upd, q2_upd, q3_upd] = updateBC(q1_upd, q2_upd, q3_upd, imn, imx, bc);
       [rho, ho, ~, a, M, p, u] = updateValues(q1_upd, q2_upd, q3_upd);
       
       [R31, R32, R33] = flux_recon(M, ho, rho, p, u, imn, imx, a, fl_rec);
       k31 = -R31/del_x; k32 = -R32/del_x; k33 = -R33/del_x;
       q1_upd(imn:imx) = q1(imn:imx) + tg*k31;
       q2_upd(imn:imx) = q2(imn:imx) + tg*k32;
       q3_upd(imn:imx) = q3(imn:imx) + tg*k33;
       [q1_upd, q2_upd, q3_upd] = updateBC(q1_upd, q2_upd, q3_upd, imn, imx, bc);
       [rho, ho, ~, a, M, p, u] = updateValues(q1_upd, q2_upd, q3_upd);
       
       [R41, R42, R43] = flux_recon(M, ho, rho, p, u, imn, imx, a, fl_rec);
       k41 = -R41/del_x; k42 = -R42/del_x; k43 = -R43/del_x;
       
       
       q1_upd(imn:imx) = q1(imn:imx) + tg*(k11 + 2*k21 + 2*k31 + k41)/6;
       q2_upd(imn:imx) = q2(imn:imx) + tg*(k12 + 2*k22 + 2*k32 + k42)/6;
       q3_upd(imn:imx) = q3(imn:imx) + tg*(k13 + 2*k23 + 2*k33 + k43)/6;
    end