function [t,pos,vel] = fload(filename,permission,datasize,datatype)
    fileid = fopen(filename,permission);
    formatSpec = datatype;
    A = fscanf(fileid,formatSpec,datasize);
    fclose(fileid);
    A = A';
    t = A(:,1);
    pos = A(:,2:4);
    vel = A(:,5:7);
end