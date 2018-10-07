function [x, del_x, imn, imx] = generate_grid(intervals)
    points =  intervals+1;
    x = linspace(-1,1,points)'; 
    del_x = x(2)-x(1);
    x = linspace(-1-del_x,1+del_x,points+2)';  % Including ghost cells
    imn = 2; imx = size(x,1)-2;        % Limits of computational domain