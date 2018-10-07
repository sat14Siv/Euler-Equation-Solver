function plotSoln(rho, u, p, T, x)
    subplot(2,2,1)
    plot(x, rho);
    ylabel('Density');
    title('Density');
    hold on;
    
    subplot(2,2,2)
    plot(x, u);
    ylabel('Velocity');
    title('Velocity');
    hold on;
    
    subplot(2,2,3)
    plot(x, p);
    ylabel('Pressure');
    xlabel('x');
    title('Pressure');
    hold on;
    
    subplot(2,2,4)
    plot(x, T);
    ylabel('Temperature');
    xlabel('x');
    title('Temperature');
    hold on;