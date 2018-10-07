function [tg, tcheck] = timeStep(CFL, u, a, tcheck, del_x)
    tinv = (max(abs(u)+ a))/del_x;      % t inverse
    tg = CFL*1/max(tinv);               % Global time step
    tcheck = tcheck+tg;