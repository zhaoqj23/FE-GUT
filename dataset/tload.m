function [t,td] = tload(filename,permission,datasize,datatype)
    fileid = fopen(filename,permission);
    formatSpec = datatype;
    A = fscanf(fileid,formatSpec,datasize);
    fclose(fileid);
    A = A';
    t = A(:,1);
    td = A(:,2);
end