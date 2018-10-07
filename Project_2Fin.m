%% Euler Equation Solver
clc; clear;

gamma = 1.4; R = 287;
%% Read Input File
filename = 'input.txt';
[CFL,fl_rec,time_integ,rhoL,rhoR,uL,uR,TL,TR,intervals,iterations,bc,...
    writeFreq] = get_input(filename);
%% Generate Grid and Initialize left and right states
[x, del_x, imn, imx] = generate_grid(intervals);
[rho, u, T, rhoWrite, uWrite, pWrite] = initialize(rhoL, rhoR, uL, uR, TL, TR,...
    bc,iterations,writeFreq, size(x,1));
%Vectors have been initilialized along with ghost cells
%% Initial Values
p = R*rho.*T;
a = (gamma*R.*T).^(0.5);
q1 = rho; q2 = rho.*u; q3 = p/(gamma - 1) + 0.5*rho.*u.^2;         
M = u./a;
ho = q3./q1 + p./rho;
%% Evolution
tcheck = 0; iter = 0; count = 0; time = 0;
while(iter < iterations ) % tcheck < tg, depending on the constraint.
    fprintf('Passed %d\n',iter);
    iter = iter + 1;

    % Flux Reconstruction
    [R1, R2, R3] = flux_recon(M, ho, rho, p, u, imn, imx, a, fl_rec);
    % Time Step Calculation
    [tg, tcheck] = timeStep(CFL, u, a, tcheck, del_x);
    % Time Integration
    [q1_upd, q2_upd, q3_upd] = timeInteg(tg, del_x, q1, q2, q3, R1, R2,...
        R3, imn, imx, time_integ,fl_rec, bc);
    % Update BCs
    [q1_upd, q2_upd, q3_upd] = updateBC(q1_upd, q2_upd, q3_upd, imn, imx, bc);
    q1 = q1_upd; q2 = q2_upd; q3 = q3_upd;
    % Sanity check
    if ((min(q1) || min(q3)) < 0)  
        fprintf('Obtaining unphysical negative values\n');
        break
    end
    % Extract primitive variables from the solution 
    [rho, ho, T, a, M, p, u] = updateValues(q1_upd, q2_upd, q3_upd);
    % Write rho, u, p
    if mod(iter, writeFreq) == 0
        [rhoWrite,uWrite,pWrite,count] = writeSolution(count,rhoWrite,...
            uWrite,pWrite,rho,u,p,iter);
    end
end
%% % Plotting
plotSoln(rho(imn:imx), u(imn:imx), p(imn:imx), T(imn:imx), x(imn+1:imx+1))

