function [CFL,fl_rec,time_integ,rhoL,rhoR,uL,uR,TL,TR,int,iters,bc,writeFreq] = get_input(file)

    fileID = fopen(file);
    C = textscan(fileID,'%f %f %f','CommentStyle','#');
    fclose(fileID);

    col1 = C{1}; col2 = C{2}; col3 = C{3};
    CFL = col1(1); fl_rec = col1(2); time_integ = col1(3); int = col1(7);
    iters = col1(8); writeFreq = col1(9); bc = col1(4);
    rhoL = col1(5); rhoR = col1(6);
    uL = col2(5); uR = col2(6);
    TL = col3(5); TR = col3(6);
    