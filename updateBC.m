function [q1_upd, q2_upd, q3_upd] = updateBC(q1_upd, q2_upd, q3_upd, imn, imx, bc)
    q1_upd(1) = q1_upd(imn);
    q1_upd(end) = q1_upd(imx);

    q2_upd(1) = bc*q2_upd(imn);
    q2_upd(end) = bc*q2_upd(imx);
    
    q3_upd(1) = q3_upd(imn);
    q3_upd(end) = q3_upd(imx);
