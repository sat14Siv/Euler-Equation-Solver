function [rhoWrite,uWrite,pWrite,count] = writeSolution(count,rhoWrite,uWrite,pWrite,rho,u,p,iter)
    count = count + 1;
    rhoWrite(2:end,count+1) = rho(2:end-1); rhoWrite(1,count+1) = iter;
    uWrite(2:end,count+1) = u(2:end-1); uWrite(1,count+1) = iter;
    pWrite(2:end,count+1) = p(2:end-1); pWrite(1,count+1) = iter;

    save('Density.txt', 'rhoWrite', '-ascii');
    save('Velocity.txt', 'uWrite', '-ascii');
    save('Pressure.txt', 'pWrite', '-ascii');